#!/usr/bin/env python
"""Damage study: what a mid-trace offset-charge jump does to the static path.

Run: ``python checks/study_charge_event_static.py``. Renders
``docs/figures/charge_event_static.png`` and prints the tables quoted in
``docs/reconstruction.md`` (charge-event section) and issue #44.

Operating point: the reference static physics scenario at n_g = 0.18
(contrast ~1.6 at 10 kSa/s, sigma 2e-4), 17 Hz background, 5 s traces, one
charge jump of amplitude delta at mid-trace. Note the parity-blind charge sits
at n_g = 0.25, only delta = 0.07 away: at this bias point a *small* jump is the
dangerous one (it lands the trace near-blind), while delta ~ 0.2-0.32 moves
each blob centre by up to ~1.3 sigma, the regime that fabricates flip clusters.
``delta = 0.5`` is the exact mirror (chi_odd(n_g) = chi_even(n_g + 1/2)): the
blob pair maps onto itself and no in-trace statistic can see the jump.

Three arms:

1. **Damage sweep** -- flip F1, fake bursts at the jump, and the decoded rate,
   with ``segment_charge_jumps`` off (before) and on (after). The after column
   is skipped with a note until the detector exists, so this script both locks
   the "before" baseline and later renders the comparison.
2. **Impact event** -- a jump coincident with an n_qp ~ 20 burst (impacts cause
   both), in-window multiplicity vs the jump-free equivalent.
3. **Null calibration** -- the maximum split-likelihood gain the detector's
   scan sees on jump-free traces; the shipped ``min_gain_nats`` must clear
   this with margin. Skipped until the detector exists.
"""
import inspect
import sys
import time
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from qpd.reconstruction import (detect_bursts,  # noqa: E402
                                reconstruct_parity_flips_static, score_flips)
from reconstruction_scenarios import build_static_scenario  # noqa: E402

NG0 = 0.18
FS = 1e4
BG = 17.0
DUR = 5.0
T_JUMP = 2.5
DELTAS = [0.02, 0.05, 0.07, 0.10, 0.15, 0.20, 0.25, 0.32, 0.40, 0.50]
N_SEEDS = 8
TOL = 0.5e-3
FIGDIR = Path(__file__).resolve().parents[1] / "docs" / "figures"
T0 = time.perf_counter()

HAVE_SEGMENTATION = "segment_charge_jumps" in inspect.signature(
    reconstruct_parity_flips_static).parameters
ARMS = [("before", {"segment_charge_jumps": False}),
        ("after", {"segment_charge_jumps": True})] if HAVE_SEGMENTATION else [
        ("before", {})]


def fake_bursts_at_jump(rec, pad=0.05):
    """Bursts overlapping the jump window; none are real (background only)."""
    bursts = detect_bursts(rec.flip_times, rec.diagnostics["rate_hz"],
                           duration=DUR)
    return sum(1 for b in bursts
               if b.t_end >= T_JUMP - pad and b.t_start <= T_JUMP + pad)


def run_damage_sweep():
    print(f"\n== 1. damage sweep: jump of delta at t = {T_JUMP} s "
          f"({len(DELTAS)} deltas x {N_SEEDS} seeds x {len(ARMS)} arm(s)) ==")
    ng_col = "n_g'"
    hdr = f"{'delta':>6} {ng_col:>5} {'contrast':>9}"
    for label, _ in ARMS:
        hdr += (f" {'F1 ' + label:>10} {'fake burst ' + label:>12} "
                f"{'rate ' + label:>10}")
    print(hdr)
    rows = {label: [] for label, _ in ARMS}
    baseline_f1 = {}
    for delta in [0.0] + DELTAS:
        scn = build_static_scenario(
            n_g=NG0, duration=DUR, sample_rate=FS, tunnel_rate_hz=BG,
            charge_jumps=(np.array([T_JUMP]), np.array([delta])))
        f1s = {label: [] for label, _ in ARMS}
        fakes = {label: 0 for label, _ in ARMS}
        rates = {label: [] for label, _ in ARMS}
        contrast = []
        for seed in range(N_SEEDS):
            tr = scn.simulate(seed=seed)
            for label, kw in ARMS:
                rec = reconstruct_parity_flips_static(tr.iq, FS, **kw)
                s = score_flips(tr.flip_times, rec.flip_times, tol=TOL)
                f1s[label].append(s.hard_f1)
                fakes[label] += fake_bursts_at_jump(rec)
                rates[label].append(rec.diagnostics["rate_hz"])
                if label == ARMS[0][0]:
                    contrast.append(rec.diagnostics["contrast"])
        line = f"{delta:6.2f} {NG0 + delta:5.2f} {np.mean(contrast):9.2f}"
        for label, _ in ARMS:
            line += (f" {np.mean(f1s[label]):10.3f} "
                     f"{fakes[label] / N_SEEDS:12.2f} "
                     f"{np.mean(rates[label]):10.1f}")
            rows[label].append((delta, float(np.mean(f1s[label])),
                                fakes[label] / N_SEEDS,
                                float(np.mean(rates[label]))))
            if delta == 0.0:
                baseline_f1[label] = float(np.mean(f1s[label]))
        print(line + ("   <- no jump (baseline)" if delta == 0.0 else ""),
              flush=True)
    print(f"({time.perf_counter() - T0:.0f}s)")
    return rows, baseline_f1


def run_impact_arm():
    print("\n== 2. impact event: n_qp ~ 20 burst at the jump time ==")
    print(f"{'delta':>6}", *(f"{'n_in_window ' + label:>18}"
                             for label, _ in ARMS))
    for delta in [0.0, 0.10, 0.20]:
        jumps = ((np.array([T_JUMP]), np.array([delta]))
                 if delta > 0 else None)
        scn = build_static_scenario(
            n_g=NG0, duration=DUR, sample_rate=FS, tunnel_rate_hz=BG,
            burst_times=np.array([T_JUMP]), burst_n_qp=20.0,
            charge_jumps=jumps)
        counts = {label: [] for label, _ in ARMS}
        n_true = []
        for seed in range(N_SEEDS):
            tr = scn.simulate(seed=100 + seed)
            if not tr.bursts:
                continue
            b = tr.bursts[0]
            lo, hi = b.t_start - 1e-3, b.t_end + 1e-3
            n_true.append(b.n_qp)
            for label, kw in ARMS:
                rec = reconstruct_parity_flips_static(tr.iq, FS, **kw)
                counts[label].append(np.count_nonzero(
                    (rec.flip_times >= lo) & (rec.flip_times <= hi)))
        line = f"{delta:6.2f}"
        for label, _ in ARMS:
            line += f" {np.mean(counts[label]):18.1f}"
        print(line + f"   (true {np.mean(n_true):.1f})", flush=True)
    print(f"({time.perf_counter() - T0:.0f}s)")


def run_null_calibration():
    print("\n== 3. null calibration: max split-gain on jump-free traces ==")
    try:
        from qpd.reconstruction.charge_events import detect_charge_events
    except ImportError:
        print("  detector not implemented yet -- skipped")
        return
    from qpd.reconstruction.static_bias import fit_two_blobs
    scn = build_static_scenario(n_g=NG0, duration=DUR, sample_rate=FS,
                                tunnel_rate_hz=BG)
    gains, n_events = [], 0
    for seed in range(50):
        tr = scn.simulate(seed=7000 + seed)
        model = fit_two_blobs(tr.iq)
        diag = detect_charge_events(tr.iq, model, 1.0 / FS)
        n_events += diag.boundaries.size
        gains.append(diag.scan.get("max_gain_nats", np.nan))
    gains = np.asarray(gains, dtype=float)
    print(f"  50 jump-free traces: {n_events} boundaries accepted, "
          f"max scan gain {np.nanmax(gains):.1f} nats "
          f"(mean {np.nanmean(gains):.1f}, p95 {np.nanpercentile(gains, 95):.1f})")
    print(f"({time.perf_counter() - T0:.0f}s)")


def render_figure(rows, baseline_f1):
    FIGDIR.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    ax = axes[0]
    for label, _ in ARMS:
        d = np.array([r[0] for r in rows[label]])
        f1 = np.array([r[1] for r in rows[label]])
        ax.plot(d[1:], f1[1:], "o-" if label == "before" else "s-",
                label=f"flip F1 ({label})")
        ax.axhline(baseline_f1[label], ls=":", lw=1,
                   color=ax.lines[-1].get_color())
    ax.axvline(0.25 - NG0, color="gray", lw=1, ls="--")
    ax.text(0.25 - NG0, 0.02, " blind landing", color="gray", fontsize=8,
            rotation=90, va="bottom")
    ax.set_xlabel(r"jump amplitude $\delta$ [$n_g$]")
    ax.set_ylabel("flip F1 (5 s trace, jump at 2.5 s)")
    ax.set_ylim(0, 1.05)
    ax.legend()
    ax.set_title(f"$n_g$ = {NG0}, {BG:.0f} Hz, 10 kSa/s (dotted: no jump)")

    ax = axes[1]
    for label, _ in ARMS:
        d = np.array([r[0] for r in rows[label]])
        fk = np.array([r[2] for r in rows[label]])
        ax.plot(d[1:], fk[1:], "o-" if label == "before" else "s-",
                label=f"fake bursts at jump ({label})")
    ax.set_xlabel(r"jump amplitude $\delta$ [$n_g$]")
    ax.set_ylabel("fake bursts per trace at the jump time")
    ax.legend()
    fig.tight_layout()
    out = FIGDIR / "charge_event_static.png"
    fig.savefig(out, dpi=150)
    print(f"\nfigure -> {out}")


def main():
    if not HAVE_SEGMENTATION:
        print("NOTE: segment_charge_jumps not implemented yet -- "
              "running the 'before' arm only.")
    rows, baseline_f1 = run_damage_sweep()
    run_impact_arm()
    run_null_calibration()
    render_figure(rows, baseline_f1)
    print(f"total {time.perf_counter() - T0:.0f}s")


if __name__ == "__main__":
    main()
