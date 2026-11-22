"""Plan, simulate, and audit LED waveforms with private physical metadata."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import time

import numpy as np
from scipy.stats import qmc

from qpd import LEDPulseModel, WhiteGaussianNoise
from qpd.mlebench.led_dataset import LEDDatasetConfig, build_led_context, _CachedLEDSimulator


PUBLIC_FIELDS = (
    "record_id", "file", "n_samples", "duration_s", "sample_rate_hz",
    "first_pulse_time_s", "on_time_s", "off_time_s",
)


def _allocate(n, weights):
    """Largest remainder, stable ties in supplied order."""
    weights = np.asarray(weights, dtype=float)
    raw = n * weights / weights.sum()
    counts = np.floor(raw).astype(int)
    order = np.argsort(-(raw - counts), kind="stable")
    counts[order[:n - counts.sum()]] += 1
    return counts


def precise_tolerance_s(sample_rate_hz):
    if not np.isfinite(sample_rate_hz) or sample_rate_hz <= 0:
        raise ValueError("sample_rate_hz must be finite and positive")
    return max(50e-6, 2.0 / sample_rate_hz)


def load_spec(path):
    spec = json.loads(Path(path).read_text())
    for key, p in spec["parameters"].items():
        lo, hi, nominal = p["min"], p["max"], p["nominal"]
        if not np.isfinite([lo, hi, nominal]).all() or not 0 <= lo <= nominal <= hi:
            raise ValueError(f"Invalid range for {key}")
        if p["sampling"] not in ("log", "linear", "log1p_positive"):
            raise ValueError(f"Unknown sampling distribution for {key}")
        if p["sampling"] == "log" and lo <= 0:
            raise ValueError(f"Nonpositive logarithmic bound for {key}")
    weights = spec["sampling_design"]["component_weights"]
    if set(weights) != {"broad", "nominal", "stress"} or any(v < 0 for v in weights.values()):
        raise ValueError("Invalid sampling components")
    if not np.isclose(sum(weights.values()), 1):
        raise ValueError("Component weights must sum to one")
    if spec["zero_led_mean_fraction"] != 0.1:
        raise ValueError("This campaign requires exactly 10% zero-LED means")
    if spec["simulation"]["veto"] or spec["scorer"]["veto"]:
        raise ValueError("No veto is allowed")
    return spec


def plan_pilot(spec, *, n_records=260, seed=11010):
    if isinstance(n_records, bool) or not isinstance(n_records, int) or n_records < 20 or n_records % 10:
        raise ValueError("Pilot size must be a multiple of ten and at least twenty")
    if isinstance(seed, bool) or not isinstance(seed, int) or seed < 0:
        raise ValueError("Seed must be a nonnegative integer")
    root = np.random.SeedSequence([seed, 731901])
    design_seed, identity_seed, realization_seed, augment_seed = root.spawn(4)
    rng, identities, realizations, augment = [np.random.default_rng(s) for s in
                                            (design_seed, identity_seed, realization_seed, augment_seed)]
    parameters = spec["parameters"]
    keys = list(parameters)
    weights = spec["sampling_design"]["component_weights"]
    components = list(weights)
    counts = _allocate(n_records, list(weights.values()))
    zero_counts = _allocate(n_records // 10, counts)
    rows = []
    # Explicit endpoint interactions, not all corners of an exponential grid.
    stress_cases = [
        ("low_rate_dense", {"sample_rate_hz": "min", "mean_tunnels": "max", "off_time_s": "min"}),
        ("noisy_resets", {"snr_max": "min", "charge_jump_rate_hz": "max", "off_time_s": "min"}),
        ("short_sparse", {"duration_s": "min", "quiescent_rate_hz": "min", "off_time_s": "max"}),
        ("maximum_record", {"sample_rate_hz": "max", "duration_s": "max"}),
        ("broad_overlap", {"burst_sigma_s": "max", "burst_tau_s": "max", "on_time_s": "max", "off_time_s": "min"}),
        ("narrow_fast_background", {"burst_sigma_s": "min", "burst_tau_s": "min", "on_time_s": "min", "quiescent_rate_hz": "max"}),
    ]
    for component, count, zero_count in zip(components, counts, zero_counts):
        # Older supported SciPy releases use seed rather than rng here.
        unit = qmc.LatinHypercube(d=len(keys), seed=int(rng.integers(2**32))).random(int(count))
        zero_indices = set(rng.choice(count, int(zero_count), replace=False).tolist())
        for index, vector in enumerate(unit):
            values = {}
            for key, u in zip(keys, vector):
                p = parameters[key]
                if p["sampling"] == "log":
                    value = np.exp(np.log(p["min"]) + u * np.log(p["max"] / p["min"]))
                elif p["sampling"] == "log1p_positive":
                    value = np.expm1(u * np.log1p(p["max"]))
                    value = max(float(value), np.nextafter(0., 1.))
                else:
                    value = p["min"] + u * (p["max"] - p["min"])
                if component == "nominal":
                    value = p["nominal"] if index % 2 == 0 else np.clip(
                        p["nominal"] * rng.uniform(.8, 1.2), p["min"], p["max"])
                values[key] = float(value)
            stress_case = None
            if component == "stress":
                stress_case, endpoints = stress_cases[index % len(stress_cases)]
                values.update({key: parameters[key][edge] for key, edge in endpoints.items()})
            if index in zero_indices:
                values["mean_tunnels"] = 0.
            values["burst_mu_s"] = spec["burst_mu_s"]
            period = values["on_time_s"] + values["off_time_s"]
            values["first_pulse_time_s"] = float(rng.uniform(0, period))
            # Ceil ensures at least the requested duration, including the 0.5s bound.
            n = int(np.ceil(values["duration_s"] * values["sample_rate_hz"]))
            n = int(np.clip(n, np.ceil(parameters["duration_s"]["min"] * values["sample_rate_hz"]),
                            np.floor(parameters["duration_s"]["max"] * values["sample_rate_hz"])))
            radius = spec["iq_augmentation"]["translation_disk_radius_in_max_contrast"] * np.sqrt(augment.random())
            angle = augment.uniform(-np.pi, np.pi)
            lo, hi = spec["iq_augmentation"]["rotation_radians"]
            rows.append(dict(record_id=identities.bytes(16).hex(), component=component,
                stress_case=stress_case, parameters=values, n_samples=n,
                duration_s=n / values["sample_rate_hz"],
                simulation_seed=int(realizations.integers(0, np.iinfo(np.int64).max)),
                iq_transform=dict(rotation_rad=float(augment.uniform(lo, hi)),
                    shift_i_in_contrast=float(radius*np.cos(angle)),
                    shift_q_in_contrast=float(radius*np.sin(angle)))))
    identities.shuffle(rows)
    assert len({row["record_id"] for row in rows}) == n_records
    assert len({row["simulation_seed"] for row in rows}) == n_records
    return rows


def transform_iq(z, transform, max_contrast):
    shift = max_contrast * complex(transform["shift_i_in_contrast"], transform["shift_q_in_contrast"])
    return np.exp(1j * transform["rotation_rad"]) * z + shift


def plan_summary(rows):
    return dict(n_records=len(rows), n_zero_led=sum(r["parameters"]["mean_tunnels"] == 0 for r in rows),
        raw_iq_bytes=8*sum(r["n_samples"] for r in rows),
        max_samples=max(r["n_samples"] for r in rows),
        total_duration_s=sum(r["duration_s"] for r in rows),
        components={c: {"records": sum(r["component"] == c for r in rows),
                          "zero_led": sum(r["component"] == c and r["parameters"]["mean_tunnels"] == 0 for r in rows)}
                    for c in ("broad", "nominal", "stress")})


def simulate_record(context, p, duration_s, seed):
    """Simulate one full trajectory with the calibrated device and specified LED settings."""
    sim = _CachedLEDSimulator(qpd=context.qpd, resonator=context.resonator,
        f_drive=context.config.f_drive_hz, sample_rate=p["sample_rate_hz"],
        gamma_even_to_odd=p["quiescent_rate_hz"], gamma_odd_to_even=p["quiescent_rate_hz"],
        noise=WhiteGaussianNoise(sigma=context.noise_sigma(p["snr_max"])),
        led_pulsing=LEDPulseModel(on_time=p["on_time_s"], off_time=p["off_time_s"],
            first_pulse_time=p["first_pulse_time_s"], expected_n_qp=p["mean_tunnels"],
            mu=p["burst_mu_s"], sigma=p["burst_sigma_s"], tau=p["burst_tau_s"],
            charge_jump_rate_hz=p["charge_jump_rate_hz"]))
    sim.cached_chi = context.chi_grid
    return sim.simulate(duration=duration_s, seed=seed)


def generate_pilot(spec, rows, destination, *, source_paths=(), scientific_split=None):
    """Write simulated records and validate their on-disk waveforms and truth."""
    if scientific_split not in (None, "train", "validation", "reserve"):
        raise ValueError("Invalid scientific split")
    purpose = "engineering_only" if scientific_split is None else "scientific_waveforms"
    out = Path(destination)
    out.mkdir(parents=True, exist_ok=False)
    (out / "INCOMPLETE").write_text("Engineering generation has not passed on-disk validation.\n")
    source_bytes = {str(p): Path(p).read_bytes() for p in source_paths}
    source_hashes = {p: hashlib.sha256(content).hexdigest() for p, content in source_bytes.items()}
    for directory in ("public/waveforms", "private/truth"):
        (out / directory).mkdir(parents=True)
    (out / "private/source_snapshot").mkdir()
    for path, content in source_bytes.items():
        (out / "private/source_snapshot" / source_hashes[path]).write_bytes(content)
    start = time.perf_counter()
    nominal = {k: v["nominal"] for k, v in spec["parameters"].items()}
    config = LEDDatasetConfig(**{k: v for k, v in nominal.items() if k not in ("snr_max", "mean_tunnels")},
        burst_mu_s=spec["burst_mu_s"], snr_max_values=(nominal["snr_max"],),
        mean_tunnels_values=(nominal["mean_tunnels"],), chunks_per_split={k: 1 for k in ("train", "validation", "test")})
    context = build_led_context(config)
    cache_seconds = time.perf_counter()-start
    write_json(out / "private/specification.json", spec)
    write_json(out / "private/design.json", rows)
    public_rows, audits = [], []
    for number, row in enumerate(rows, 1):
        begin = time.perf_counter()
        p = row["parameters"]
        result = simulate_record(context, p, row["duration_s"], row["simulation_seed"])
        z = transform_iq(result.iq, row["iq_transform"], context.max_separation)
        i, q = z.real.astype(np.float32), z.imag.astype(np.float32)
        if len(i) != row["n_samples"] or not np.isfinite(i).all() or not np.isfinite(q).all():
            raise AssertionError("Invalid waveform")
        truth, led = result.flip_times, result.led
        if not np.all((truth >= 0) & (truth <= result.t[-1])) or np.any(np.diff(truth) < 0):
            raise AssertionError("Invalid flip truth")
        events = led.event_times
        if p["mean_tunnels"] == 0 and len(events):
            raise AssertionError("Zero-LED record contains LED-induced flips")
        for burst in led.bursts:
            if np.any(burst.event_times < burst.t_arrival):
                raise AssertionError("Noncausal LED event")
        file = f"waveforms/{row['record_id']}.npz"
        np.savez_compressed(out / "public" / file, I=i, Q=q)
        np.savez_compressed(out / "private/truth" / f"{row['record_id']}.npz",
            flip_times=truth, pulse_times=led.pulse_times, pulse_end_times=led.pulse_end_times,
            led_event_times=events, charge_reset_times=led.charge_reset_times,
            charge_reset_values=led.charge_reset_values, charge_reset_is_led=led.charge_reset_is_led,
            initial_n_g=led.initial_n_g)
        public_rows.append(dict(record_id=row["record_id"], file=file, n_samples=len(i),
            duration_s=row["duration_s"], **{key: p[key] for key in
                ("sample_rate_hz", "first_pulse_time_s", "on_time_s", "off_time_s")}))
        bins = np.ceil(truth * p["sample_rate_hz"]).astype(np.int64)
        _, interval_counts = np.unique(bins, return_counts=True)
        audits.append(dict(record_id=row["record_id"], n_truth=len(truth),
            n_led_pulses=len(led.pulse_times), n_led_events_full_pulses=len(events),
            n_led_events_recorded=int(np.count_nonzero(events <= result.t[-1])),
            n_additional_charge_resets=int(np.count_nonzero(~led.charge_reset_is_led)),
            n_truth_in_multiply_occupied_sample_intervals=int(interval_counts[interval_counts > 1].sum()),
            precise_tolerance_s=precise_tolerance_s(p["sample_rate_hz"]),
            last_sample_s=float(result.t[-1]),
            waveform_sha256=hashlib.sha256(i.tobytes()+q.tobytes()).hexdigest(),
            generation_and_write_seconds=time.perf_counter()-begin))
        if number % 20 == 0 or number == len(rows):
            print(json.dumps(dict(stage="pilot", completed=number, total=len(rows))), flush=True)
    hashes = [r["waveform_sha256"] for r in audits]
    if len(set(hashes)) != len(hashes):
        raise AssertionError("Repeated waveform bytes")
    write_json(out / "public/manifest.json", public_rows)
    write_json(out / "public/dataset.json", dict(version=spec["version"], purpose=purpose,
        training_ready=False, baseline_features="not_computed", veto=False,
        observation_window="first through last sample", model_metadata=list(PUBLIC_FIELDS),
        nominal_emg=spec["baseline"]["nominal_emg"]))
    summary = dict(**plan_summary(rows), cache_seconds=cache_seconds,
        elapsed_seconds=time.perf_counter()-start, max_parity_contrast=context.max_separation,
        compressed_iq_bytes=sum(p.stat().st_size for p in (out / "public/waveforms").glob("*.npz")),
        n_truth=sum(r["n_truth"] for r in audits),
        n_truth_in_multiply_occupied_sample_intervals=sum(r["n_truth_in_multiply_occupied_sample_intervals"] for r in audits),
        training_ready=False, model_or_hmm_execution=False,
        generated_splits=[scientific_split or "engineering"], source_sha256=source_hashes)
    write_json(out / "private/audit.json", audits)
    write_json(out / "private/summary.json", summary)
    audit_pilot(out)
    return summary


def audit_pilot(destination):
    """Validate every exported file before issuing an engineering completion receipt."""
    out = Path(destination)
    (out / "COMPLETE.json").unlink(missing_ok=True)
    (out / "INCOMPLETE").write_text("On-disk engineering validation is in progress or failed.\n")
    rows = json.loads((out / "public/manifest.json").read_text())
    design = json.loads((out / "private/design.json").read_text())
    audits = json.loads((out / "private/audit.json").read_text())
    summary = json.loads((out / "private/summary.json").read_text())
    metadata = json.loads((out / "public/dataset.json").read_text())
    if metadata["training_ready"] or metadata["purpose"] not in ("engineering_only", "scientific_waveforms") or metadata["veto"]:
        raise ValueError("Invalid engineering status")
    expected = {r["record_id"]: r for r in design}
    hashes = {r["record_id"]: r["waveform_sha256"] for r in audits}
    ids = [r["record_id"] for r in rows]
    if len(set(ids)) != len(ids) or set(ids) != set(expected) or set(ids) != set(hashes):
        raise ValueError("Record identities differ")
    if len(rows) != summary["n_records"] or sum(r["parameters"]["mean_tunnels"] == 0 for r in design)*10 != len(rows):
        raise ValueError("Incorrect record allocation")
    permitted = {"manifest.json", "dataset.json"} | {f"waveforms/{i}.npz" for i in ids}
    present = {p.relative_to(out / "public").as_posix() for p in (out / "public").rglob("*") if p.is_file()}
    if present != permitted:
        raise ValueError("Unexpected or missing public artifact")
    if {p.name for p in (out / "private/truth").iterdir()} != {f"{i}.npz" for i in ids}:
        raise ValueError("Unexpected or missing truth artifact")
    for row in rows:
        if set(row) != set(PUBLIC_FIELDS) or row["file"] != f"waveforms/{row['record_id']}.npz":
            raise ValueError("Public schema or path violation")
        d = expected[row["record_id"]]
        for key in ("sample_rate_hz", "first_pulse_time_s", "on_time_s", "off_time_s"):
            if row[key] != d["parameters"][key]:
                raise ValueError("Acquisition metadata mismatch")
        with np.load(out / "public" / row["file"], allow_pickle=False) as data:
            if set(data.files) != {"I", "Q"}:
                raise ValueError("Unexpected waveform channel")
            i, q = data["I"], data["Q"]
            if i.dtype != np.float32 or q.dtype != np.float32 or i.shape != (row["n_samples"],) or q.shape != i.shape:
                raise ValueError("Waveform shape or dtype mismatch")
            if not np.isfinite(i).all() or not np.isfinite(q).all():
                raise ValueError("Nonfinite waveform")
            if hashlib.sha256(i.tobytes()+q.tobytes()).hexdigest() != hashes[row["record_id"]]:
                raise ValueError("Waveform hash mismatch")
        if row["n_samples"] != d["n_samples"] or row["duration_s"] != d["duration_s"]:
            raise ValueError("Recording duration mismatch")
        with np.load(out / "private/truth" / f"{row['record_id']}.npz", allow_pickle=False) as data:
            times = data["flip_times"]
            if not np.isfinite(times).all() or np.any(np.diff(times) < 0) or np.any(times < 0) or np.any(times > (row["n_samples"]-1)/row["sample_rate_hz"]):
                raise ValueError("Truth outside observed window")
            if d["parameters"]["mean_tunnels"] == 0 and data["led_event_times"].size:
                raise ValueError("Zero LED intensity emitted LED events")
    receipt = dict(status="complete_" + metadata["purpose"], n_records=len(rows), n_zero_led=len(rows)//10,
        checks=["public allowlist", "identities and allocation", "all waveform shapes/dtypes/hashes",
                "observable metadata", "truth boundaries", "zero LED semantics"],
        manifest_sha256=hashlib.sha256((out / "public/manifest.json").read_bytes()).hexdigest(),
        training_ready=False)
    write_json(out / "private/validation.json", receipt)
    write_json(out / "COMPLETE.json", receipt)
    (out / "INCOMPLETE").unlink(missing_ok=True)
    return receipt


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
