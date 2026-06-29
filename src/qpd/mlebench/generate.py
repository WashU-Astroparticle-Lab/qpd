"""Generate MLE-bench-style QPD datasets of complex S21 I/Q time traces.

Each *chunk* is an independent, self-contained ``chunk_seconds``-long trace
sampled at ``physics.sample_rate_hz``. Chunks share one frozen
:class:`~qpd.mlebench.config.PhysicsConfig`; per chunk we randomise only the
*structural* events (sawtooth phase, charge-jump times/sizes, burst onsets) and
the simulator's noise/telegraph seed. The signal dataset enables bursts; the
background dataset disables them. Everything else is identical.

Layout written per competition (``qpd-signal`` / ``qpd-background``)::

    <competition_id>/
        description.md
        metadata.json
        manifest.csv
        prepared/
            public/
                train/chunk_XXXX.parquet      # I, Q (float32)
                test/chunk_XXXX.parquet       # I, Q (float32)
                train_labels.csv              # long-format truth for TRAIN
                sample_submission.csv         # submission schema
            private/
                test_labels.csv               # long-format truth for TEST (key)

The time axis is uniform and reconstructed from the manifest:
``t = t0 + arange(n_samples) / sample_rate`` (``t0 = 0`` per chunk).
"""

from __future__ import annotations

import json
import math
import os
import time
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from qpd import (
    QPD,
    ChargeJumpEvents,
    QuasiparticleBurstModel,
    ResonatorConfig,
    SawtoothNg,
    VNASimulator,
    WhiteGaussianNoise,
    poisson_burst_times,
)

from qpd.mlebench.config import (
    DEFAULT_PHYSICS,
    DatasetConfig,
    GenerationBudget,
    PhysicsConfig,
    default_dataset_configs,
)

# Long-format event-table schema shared by labels, answers and submissions.
LABEL_COLUMNS = ["chunk_id", "event_type", "t", "value"]
_DATASET_TAG = {"signal": 0, "background": 1}


@dataclass
class ChunkTruth:
    """Ground truth for a single chunk."""

    chunk_id: str
    dataset: str
    duration: float
    sample_rate: float
    n_samples: int
    seed: int
    # Parity-flip (quasiparticle tunneling) times -- background telegraph and,
    # for the signal dataset, burst events merged. Exact, sub-sample.
    flip_times: np.ndarray
    # Charge jumps: time and fractional jump size (wrapped to (-0.5, 0.5]).
    charge_times: np.ndarray
    charge_fracs: np.ndarray
    # Bursts (signal only): onset time (EMG reference), QP count, event span.
    burst_onsets: np.ndarray = field(default_factory=lambda: np.empty(0))
    burst_n_qp: np.ndarray = field(default_factory=lambda: np.empty(0, int))
    burst_t_start: np.ndarray = field(default_factory=lambda: np.empty(0))
    burst_t_end: np.ndarray = field(default_factory=lambda: np.empty(0))

    def label_rows(self) -> list[dict]:
        rows: list[dict] = []
        for t in self.flip_times:
            rows.append(
                {"chunk_id": self.chunk_id, "event_type": "flip",
                 "t": float(t), "value": math.nan}
            )
        for t, v in zip(self.charge_times, self.charge_fracs):
            rows.append(
                {"chunk_id": self.chunk_id, "event_type": "charge",
                 "t": float(t), "value": float(v)}
            )
        for t, n in zip(self.burst_onsets, self.burst_n_qp):
            rows.append(
                {"chunk_id": self.chunk_id, "event_type": "burst",
                 "t": float(t), "value": float(n)}
            )
        return rows


# --------------------------------------------------------------------------- #
# Simulator context (built once, reused across chunks)
# --------------------------------------------------------------------------- #


@dataclass
class _SimContext:
    qpd: QPD
    resonator: ResonatorConfig
    f_drive: float
    physics: PhysicsConfig


def build_context(physics: PhysicsConfig = DEFAULT_PHYSICS) -> _SimContext:
    """Build the QPD device + resonator once (bare-cavity back-solve is reused)."""
    qpd = QPD(e_j_hz=physics.e_j_hz, e_c_hz=physics.e_c_hz)
    qpd.coupling_g_hz = physics.coupling_g_hz
    f_bare = physics.solve_bare_fr_hz(qpd)
    resonator = ResonatorConfig(
        f_r=f_bare,
        q_i=physics.q_i,
        q_c_abs=physics.q_c_abs,
        phi=physics.phi,
        a=physics.a,
        alpha=physics.alpha,
        tau=physics.tau_env,
    )
    return _SimContext(
        qpd=qpd,
        resonator=resonator,
        f_drive=physics.f_drive_hz(),
        physics=physics,
    )


def _wrap_frac(x: np.ndarray) -> np.ndarray:
    """Wrap charge fractions to (-0.5, 0.5]."""
    x = np.asarray(x, dtype=float)
    return ((x + 0.5) % 1.0) - 0.5


# --------------------------------------------------------------------------- #
# Single chunk
# --------------------------------------------------------------------------- #


def generate_chunk(
    dcfg: DatasetConfig,
    chunk_index: int,
    master_seed: int,
    context: _SimContext | None = None,
) -> tuple[np.ndarray, np.ndarray, ChunkTruth]:
    """Generate one chunk: return (I_float32, Q_float32, truth).

    Deterministic in ``(master_seed, dataset, chunk_index)``.
    """
    if context is None:
        context = build_context()
    physics = context.physics
    duration = float(dcfg.chunk_seconds)
    tag = _DATASET_TAG.get(dcfg.name, abs(hash(dcfg.name)) % 1000)

    # Independent RNG streams: structural draws vs simulator (noise/telegraph).
    struct_rng = np.random.default_rng(
        np.random.SeedSequence([master_seed, tag, chunk_index, 1])
    )
    sim_seed = int(
        np.random.SeedSequence(
            [master_seed, tag, chunk_index, 2]
        ).generate_state(1)[0]
    )

    # --- Slow sawtooth drift with a random phase ---
    span = physics.sawtooth_n_g_max - physics.sawtooth_n_g_min
    period = span / abs(physics.sawtooth_slope_per_s)
    t0 = float(struct_rng.uniform(0.0, period))
    offset = SawtoothNg(
        n_g_min=physics.sawtooth_n_g_min,
        n_g_max=physics.sawtooth_n_g_max,
        slope=physics.sawtooth_slope_per_s,
        t0=t0,
    )

    # --- Charge-jump events (fractional jumps only; integer part unobservable)
    n_jumps = int(struct_rng.poisson(physics.charge_jump_rate_hz * duration))
    if n_jumps > 0:
        cj_times = np.sort(struct_rng.uniform(0.0, duration, n_jumps))
        cj_fracs = struct_rng.uniform(
            physics.charge_jump_frac_low,
            physics.charge_jump_frac_high,
            n_jumps,
        )
        offset = offset + ChargeJumpEvents(times=cj_times, deltas=cj_fracs)
    else:
        cj_times = np.empty(0, dtype=float)
        cj_fracs = np.empty(0, dtype=float)

    # --- Quasiparticle bursts (signal dataset only) ---
    qp_bursts = None
    if dcfg.bursts_enabled and physics.burst_onset_rate_hz > 0:
        onsets = poisson_burst_times(
            physics.burst_onset_rate_hz, duration, struct_rng
        )
        if onsets.size:
            qp_bursts = QuasiparticleBurstModel(
                times=onsets,
                tau=physics.burst_tau,
                mu=physics.burst_mu,
                sigma=physics.burst_sigma,
                expected_n_qp=physics.burst_expected_n_qp,
            )

    sim = VNASimulator(
        qpd=context.qpd,
        resonator=context.resonator,
        f_drive=context.f_drive,
        sample_rate=physics.sample_rate_hz,
        gamma_even_to_odd=physics.gamma_even_to_odd_hz,
        gamma_odd_to_even=physics.gamma_odd_to_even_hz,
        noise=WhiteGaussianNoise(sigma=physics.noise_sigma),
        offset_charge=offset,
        quasiparticle_bursts=qp_bursts,
        n_g_grid_points=physics.n_g_grid_points,
        charge_cutoff=physics.charge_cutoff,
    )
    result = sim.simulate(duration=duration, seed=sim_seed)

    # Burst truth (only kept inside the chunk window).
    onsets = np.array([b.t_arrival for b in result.bursts], dtype=float)
    n_qp = np.array([b.n_qp for b in result.bursts], dtype=int)
    bstart = np.array([b.t_start for b in result.bursts], dtype=float)
    bend = np.array([b.t_end for b in result.bursts], dtype=float)
    keep = (onsets >= 0.0) & (onsets < duration) if onsets.size else slice(0, 0)

    truth = ChunkTruth(
        chunk_id=f"{dcfg.name}_{chunk_index:04d}",
        dataset=dcfg.name,
        duration=duration,
        sample_rate=physics.sample_rate_hz,
        n_samples=int(result.t.size),
        seed=sim_seed,
        flip_times=np.sort(result.flip_times.astype(float)),
        charge_times=cj_times,
        charge_fracs=_wrap_frac(cj_fracs),
        burst_onsets=onsets[keep],
        burst_n_qp=n_qp[keep],
        burst_t_start=bstart[keep],
        burst_t_end=bend[keep],
    )
    return result.i.astype(np.float32), result.q.astype(np.float32), truth


# --------------------------------------------------------------------------- #
# Budget planning
# --------------------------------------------------------------------------- #


def estimate_bytes_per_chunk(
    chunk_seconds: float, sample_rate_hz: float, overhead: float = 1.02
) -> float:
    """Estimate Parquet bytes/chunk: I+Q float32 (4 B each) * n_samples."""
    n = int(round(chunk_seconds * sample_rate_hz))
    return n * 4 * 2 * overhead


def plan_generation(
    budget: GenerationBudget,
    chunk_seconds: float,
    sample_rate_hz: float,
    n_datasets: int = 2,
    sec_per_chunk: float = 1.0,
    bytes_per_chunk: float | None = None,
) -> dict:
    """Compute the per-dataset chunk count so both datasets are equal length.

    Each dataset gets ``1 / n_datasets`` of the total byte and compute budgets;
    the binding constraint wins.
    """
    if bytes_per_chunk is None:
        bytes_per_chunk = estimate_bytes_per_chunk(chunk_seconds, sample_rate_hz)
    byte_cap = (budget.max_total_bytes / n_datasets) / bytes_per_chunk
    compute_cap = (
        budget.max_total_compute_seconds / n_datasets
    ) / sec_per_chunk
    n = math.floor(min(byte_cap, compute_cap))
    if budget.max_chunks_per_dataset is not None:
        n = min(n, budget.max_chunks_per_dataset)
    n = max(n, 0)
    binding = "bytes" if byte_cap <= compute_cap else "compute"
    if budget.max_chunks_per_dataset is not None and n == budget.max_chunks_per_dataset:
        binding = "max_chunks"
    return {
        "chunks_per_dataset": int(n),
        "n_datasets": n_datasets,
        "bytes_per_chunk": float(bytes_per_chunk),
        "sec_per_chunk": float(sec_per_chunk),
        "binding_constraint": binding,
        "projected_total_bytes": float(n * n_datasets * bytes_per_chunk),
        "projected_total_seconds": float(n * n_datasets * sec_per_chunk),
        "projected_data_seconds_total": float(n * n_datasets * chunk_seconds),
    }


# --------------------------------------------------------------------------- #
# Full dataset
# --------------------------------------------------------------------------- #


def _write_chunk_parquet(path: str, i: np.ndarray, q: np.ndarray) -> int:
    df = pd.DataFrame({"I": i, "Q": q})
    df.to_parquet(path, engine="pyarrow", compression="snappy", index=False)
    return os.path.getsize(path)


def generate_dataset(
    out_root: str,
    dcfg: DatasetConfig,
    n_chunks: int,
    master_seed: int,
    physics: PhysicsConfig = DEFAULT_PHYSICS,
    context: _SimContext | None = None,
    byte_budget: float | None = None,
    progress: bool = True,
) -> dict:
    """Generate one competition dataset. Returns a summary dict."""
    if context is None:
        context = build_context(physics)

    comp_dir = os.path.join(out_root, dcfg.competition_id)
    public = os.path.join(comp_dir, "prepared", "public")
    private = os.path.join(comp_dir, "prepared", "private")
    train_dir = os.path.join(public, "train")
    test_dir = os.path.join(public, "test")
    for d in (train_dir, test_dir, private):
        os.makedirs(d, exist_ok=True)

    n_train = int(round(dcfg.train_fraction * n_chunks))
    train_label_rows: list[dict] = []
    test_label_rows: list[dict] = []
    manifest_rows: list[dict] = []

    bytes_written = 0
    t_start = time.time()
    produced = 0
    for idx in range(n_chunks):
        i, q, truth = generate_chunk(dcfg, idx, master_seed, context=context)
        split = "train" if idx < n_train else "test"
        dest = train_dir if split == "train" else test_dir
        fname = f"chunk_{idx:04d}.parquet"
        size = _write_chunk_parquet(os.path.join(dest, fname), i, q)
        bytes_written += size

        rows = truth.label_rows()
        (train_label_rows if split == "train" else test_label_rows).extend(rows)
        manifest_rows.append(
            {
                "chunk_id": truth.chunk_id,
                "split": split,
                "file": os.path.join("prepared", "public", split, fname),
                "n_samples": truth.n_samples,
                "sample_rate_hz": truth.sample_rate,
                "t0": 0.0,
                "duration_s": truth.duration,
                "seed": truth.seed,
                "n_flips": int(truth.flip_times.size),
                "n_charge": int(truth.charge_times.size),
                "n_bursts": int(truth.burst_onsets.size),
                "bytes": size,
            }
        )
        produced += 1
        if progress:
            print(
                f"  [{dcfg.name}] chunk {idx + 1}/{n_chunks} "
                f"({split}) {size / 1e6:.1f} MB  "
                f"cum {bytes_written / 1e6:.0f} MB",
                flush=True,
            )
        if byte_budget is not None and bytes_written >= byte_budget:
            if progress:
                print(
                    f"  [{dcfg.name}] byte budget reached at chunk "
                    f"{idx + 1}; stopping early.",
                    flush=True,
                )
            break

    elapsed = time.time() - t_start

    # --- Write label / submission / manifest / metadata ---
    empty = pd.DataFrame(columns=LABEL_COLUMNS)
    pd.DataFrame(train_label_rows or None, columns=LABEL_COLUMNS).to_csv(
        os.path.join(public, "train_labels.csv"), index=False
    )
    pd.DataFrame(test_label_rows or None, columns=LABEL_COLUMNS).to_csv(
        os.path.join(private, "test_labels.csv"), index=False
    )
    empty.to_csv(os.path.join(public, "sample_submission.csv"), index=False)
    pd.DataFrame(manifest_rows).to_csv(
        os.path.join(comp_dir, "manifest.csv"), index=False
    )

    summary = {
        "competition_id": dcfg.competition_id,
        "dataset": dcfg.name,
        "bursts_enabled": dcfg.bursts_enabled,
        "chunks_produced": produced,
        "n_train": min(n_train, produced),
        "n_test": max(produced - n_train, 0),
        "chunk_seconds": dcfg.chunk_seconds,
        "bytes_written": bytes_written,
        "elapsed_seconds": elapsed,
        "master_seed": master_seed,
    }
    metadata = {
        "summary": summary,
        "physics": physics.to_dict(),
        "label_schema": LABEL_COLUMNS,
        "time_axis": "t = t0 + arange(n_samples) / sample_rate_hz (t0 = 0)",
    }
    with open(os.path.join(comp_dir, "metadata.json"), "w") as f:
        json.dump(metadata, f, indent=2)
    _write_description(comp_dir, dcfg, physics, summary)
    return summary


def generate_all(
    out_root: str,
    budget: GenerationBudget = GenerationBudget(),
    physics: PhysicsConfig = DEFAULT_PHYSICS,
    chunk_seconds: float = 30.0,
    train_fraction: float = 0.6,
    master_seed: int = 0,
    sec_per_chunk: float = 1.0,
    bytes_per_chunk: float | None = None,
    progress: bool = True,
) -> dict:
    """Generate BOTH datasets (signal + background) at equal chunk count."""
    context = build_context(physics)
    plan = plan_generation(
        budget,
        chunk_seconds=chunk_seconds,
        sample_rate_hz=physics.sample_rate_hz,
        n_datasets=2,
        sec_per_chunk=sec_per_chunk,
        bytes_per_chunk=bytes_per_chunk,
    )
    n = plan["chunks_per_dataset"]
    if progress:
        print(f"Plan: {n} chunks/dataset, binding = {plan['binding_constraint']}, "
              f"projected {plan['projected_total_bytes'] / 1e9:.2f} GB, "
              f"{plan['projected_data_seconds_total'] / 60:.1f} min of data.",
              flush=True)
    per_dataset_byte_budget = budget.max_total_bytes / 2
    results = []
    for dcfg in default_dataset_configs(chunk_seconds, train_fraction):
        results.append(
            generate_dataset(
                out_root,
                dcfg,
                n_chunks=n,
                master_seed=master_seed,
                physics=physics,
                context=context,
                byte_budget=per_dataset_byte_budget,
                progress=progress,
            )
        )
    return {"plan": plan, "datasets": results}


# --------------------------------------------------------------------------- #
# description.md
# --------------------------------------------------------------------------- #


def _write_description(comp_dir, dcfg, physics, summary) -> None:
    burst_line = (
        "**Quasiparticle bursts ARE present** in this dataset "
        f"(onset rate ~{physics.burst_onset_rate_hz} Hz, "
        f"Poisson({physics.burst_expected_n_qp:.0f}) QP per burst, "
        "EMG arrival profile)."
        if dcfg.bursts_enabled
        else "**No quasiparticle bursts** are present in this dataset "
        "(background condition only)."
    )
    text = f"""# {dcfg.competition_id}: QPD parity / charge / burst reconstruction

## Overview

Complex transmission **S21 = I + iQ** time traces from a simulated **QPD
transmon** dispersively coupled to a notch readout resonator. Each example is a
self-contained **{dcfg.chunk_seconds:.0f} s** chunk sampled at
**{physics.sample_rate_hz:.0f} Hz** (uniform grid).

Both the `qpd-signal` and `qpd-background` competitions share the same
background condition: electronic white noise, a quiescent two-state
quasiparticle tunneling telegraph, and discrete charge-jump events.
{burst_line}

## Task

For each test chunk, reconstruct the ground-truth event list:

1. **Parity flips** — the time of every quasiparticle tunneling event (each
   flips the charge parity even↔odd).
2. **Charge jumps** — the time and the **fractional part** of each offset-charge
   jump (wrapped to (-0.5, 0.5]; the integer part is physically unobservable).
3. **Quasiparticle bursts** — the **onset time** and the **number of QP
   tunneling events** of each burst (this dataset only contains bursts if it is
   `qpd-signal`).

## Data

```
prepared/public/train/chunk_XXXX.parquet   # I, Q  (float32)
prepared/public/test/chunk_XXXX.parquet    # I, Q  (float32)  -- predict these
prepared/public/train_labels.csv           # ground truth for the train chunks
prepared/public/sample_submission.csv      # submission schema
manifest.csv                               # per-chunk sample_rate, n_samples, t0
```

Each Parquet file has two columns, `I` and `Q` (float32). The time stamps are a
uniform grid reconstructed from `manifest.csv`:

```
t = t0 + arange(n_samples) / sample_rate_hz      # t0 = 0
```

## Submission format

A single long-format CSV with one row per predicted event and columns:

| column     | meaning |
|------------|---------|
| `chunk_id` | e.g. `{dcfg.name}_0042` |
| `event_type` | `flip`, `charge`, or `burst` |
| `t`        | event time within the chunk [s] |
| `value`    | `flip`: blank; `charge`: fractional jump; `burst`: QP count |

Predict no rows for a chunk to predict "no events". See
`prepared/public/sample_submission.csv`.

## Grading

Per chunk and per event type, predicted events are matched to truth events
**optimally** (minimum total timing error, Hungarian assignment) within a
per-type time tolerance. Each matched event then earns continuous partial
credit from its timing residual, `q = exp(-(dt / scale)^2)`, and a *soft* F1
folds this into detection: `soft_TP = sum(q)`, `soft_precision = soft_TP /
n_pred`, `soft_recall = soft_TP / n_truth`. So both false positives/negatives
and loose timing reduce the score continuously.

- `flip`  : soft timing-F1 (tolerance 0.5 ms).
- `charge`: soft timing-F1 (tolerance 5 ms) + circular MAE of the fractional value.
- `burst` : soft timing-F1 (tolerance 5 ms) + MAE of the QP count.

The overall score is the mean soft timing-F1 across the event types present,
lightly penalised by the normalised value error. See `qpd.mlebench.grade`
(`grade` for the scalar, `grade_report(..., return_matches=True)` for the full
per-event breakdown). Higher is better.

## Provenance

Generated by `qpd.mlebench` from `notebooks/simulation.ipynb` physics.
Chunks: {summary['chunks_produced']} ({summary['n_train']} train /
{summary['n_test']} test). Master seed: {summary['master_seed']}.
Full physics parameters are in `metadata.json`.
"""
    with open(os.path.join(comp_dir, "description.md"), "w") as f:
        f.write(text)
