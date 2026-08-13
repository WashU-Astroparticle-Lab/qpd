#!/usr/bin/env python
"""Before/after study: two-state HMM vs burst-aware parity x regime HMM.

Run: ``python checks/study_burst_aware_comparison.py``. Renders
``docs/figures/burst_aware_comparison.png`` and prints the comparison
tables quoted in ``docs/reconstruction.md`` section 12f and issue #40.

Operating point (matches the issue): 10 kSa/s, contrast 2.4, 17 Hz
background, EMG burst profile tau = 1 ms / mu = 1.2 ms / sigma = 0.4 ms,
``detect_bursts`` defaults (3 ms linking, min 3 flips, no p-value gate).
"""
import sys
import time
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from qpd.reconstruction import (detect_bursts,  # noqa: E402
                                reconstruct_parity_flips_static, score_flips,
                                sweep_burst_size)
from qpd.reconstruction.benchmark import TraceFidelity  # noqa: E402
from qpd.reconstruction.static_bias import StaticBlobModel  # noqa: E402
from qpd.simulator import QuasiparticleBurstModel  # noqa: E402

FS = 1e4
CONTRAST = 2.4
BG = 17.0
FIGDIR = Path(__file__).resolve().parents[1] / "docs" / "figures"
T0 = time.perf_counter()


def make_fidelity(**recon_kwargs):
    blob = StaticBlobModel(direction=1 + 0j, origin=0j,
                           mu_a=+CONTRAST / 2, mu_b=-CONTRAST / 2,
                           sigma=1.0, weight_a=0.5)
    return TraceFidelity(
        mode="static", sample_rate=FS, n_samples=int(FS), t0=0.0, sigma=1.0,
        contrast_median=CONTRAST, contrast_max=CONTRAST, rate_hz=BG,
        n_flips=0, degenerate=False, sample_fidelity=0.9,
        decoded_fidelity=0.99, model=blob, recon_kwargs=recon_kwargs)


def main():
    fid_base = make_fidelity()
    fid_ba = make_fidelity(burst_aware=True)

    # -- burst-size sweep, both decoders --------------------------------
    nqp = [2, 3, 5, 8, 12, 20, 30, 50, 80]
    print("sweep_burst_size (before) ...", flush=True)
    pts_base = sweep_burst_size(fid_base, nqp, n_trials=30, seed=100)
    print(f"sweep_burst_size (after) ... ({time.perf_counter()-T0:.0f}s)",
          flush=True)
    pts_ba = sweep_burst_size(fid_ba, nqp, n_trials=30, seed=100)

    print(f"\n{'n_qp':>5} {'eff before':>12} {'eff after':>12} "
          f"{'<det> before':>13} {'<det> after':>12}")
    for pb, pa in zip(pts_base, pts_ba):
        print(f"{pb.n_qp_expected:5.0f} "
              f"{pb.efficiency:7.3f}+-{pb.efficiency_err:.3f} "
              f"{pa.efficiency:7.3f}+-{pa.efficiency_err:.3f} "
              f"{pb.mean_n_qp_detected:13.1f} {pa.mean_n_qp_detected:12.1f}")

    # -- background guard -----------------------------------------------
    print(f"\nbackground guard ... ({time.perf_counter()-T0:.0f}s)",
          flush=True)
    for label, kw in [("before", {"burst_aware": False}),
                      ("after", {"burst_aware": True})]:
        n_false, effs, purs = 0, [], []
        for k in range(40):
            iq, truth = fid_base.synthesize(seed=5000 + k, rate_hz=BG)
            rec = reconstruct_parity_flips_static(iq, FS, **kw)
            n_false += len(detect_bursts(rec.flip_times, BG, duration=1.0))
            s = score_flips(truth, rec.flip_times, tol=0.5e-3)
            effs.append(s.efficiency)
            purs.append(s.purity)
        print(f"  {label:6s} false bursts/s = {n_false / 40:.3f}  "
              f"flip eff = {np.nanmean(effs):.3f}  "
              f"purity = {np.nanmean(purs):.3f}")

    # -- single-burst window recovery vs the visible ceiling ------------
    print(f"\nwindow recovery ... ({time.perf_counter()-T0:.0f}s)",
          flush=True)
    table = []
    for target in [8, 20, 80]:
        model = QuasiparticleBurstModel(times=np.array([0.5]), tau=1e-3,
                                        mu=1.2e-3, sigma=0.4e-3,
                                        expected_n_qp=float(target))
        vis, dec_b, dec_a, tru = [], [], [], []
        for k in range(30):
            iq, flips, bt = fid_base.synthesize(
                seed=9000 + 31 * target + k, rate_hz=BG, bursts=model,
                return_bursts=True)
            b = bt[0]
            if b.n_qp <= 0:
                continue
            lo, hi = b.t_start - 1e-3, b.t_end + 1e-3
            t = np.arange(int(FS)) / FS
            par = np.searchsorted(flips, t, side="right") % 2
            w = (t >= lo) & (t <= hi)
            vis.append(np.count_nonzero(np.diff(par[w.nonzero()[0]])))
            tru.append(b.n_qp)
            for kw, out in [({"burst_aware": False}, dec_b),
                           ({"burst_aware": True}, dec_a)]:
                rec = reconstruct_parity_flips_static(iq, FS, **kw)
                out.append(np.count_nonzero((rec.flip_times >= lo)
                                            & (rec.flip_times <= hi)))
        table.append((float(np.mean(tru)), float(np.mean(vis)),
                      float(np.mean(dec_b)), float(np.mean(dec_a))))
        print(f"  n_qp~{target:3d}: true {table[-1][0]:5.1f}  visible "
              f"{table[-1][1]:5.1f}  before {table[-1][2]:5.1f}  "
              f"after {table[-1][3]:5.1f}")

    # -- figure ----------------------------------------------------------
    FIGDIR.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    ax = axes[0]
    x = [p.n_qp_expected for p in pts_base]
    ax.errorbar(x, [p.efficiency for p in pts_base],
                yerr=[p.efficiency_err for p in pts_base],
                marker="o", capsize=3, label="two-state HMM (before)")
    ax.errorbar(x, [p.efficiency for p in pts_ba],
                yerr=[p.efficiency_err for p in pts_ba],
                marker="s", capsize=3, label="burst-aware HMM (after)")
    ax.set_xscale("log")
    ax.set_xlabel("expected burst quasiparticle number")
    ax.set_ylabel("burst detection efficiency")
    ax.set_ylim(0, 1.05)
    ax.legend()
    ax.set_title(f"contrast {CONTRAST}, {BG:.0f} Hz background, 10 kSa/s")

    ax = axes[1]
    ax.plot([1, 100], [1, 100], "k:", lw=1, label="detected = true")
    ax.plot([p.mean_n_qp_true for p in pts_base],
            [p.mean_n_qp_detected for p in pts_base],
            "o-", label="two-state HMM (before)")
    ax.plot([p.mean_n_qp_true for p in pts_ba],
            [p.mean_n_qp_detected for p in pts_ba],
            "s-", label="burst-aware HMM (after)")
    ax.plot([r[0] for r in table], [r[1] for r in table], "^--", color="gray",
            label="visible on 10 kSa/s grid (ceiling)")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("true burst quasiparticle number")
    ax.set_ylabel("mean detected multiplicity")
    ax.legend()
    fig.tight_layout()
    out = FIGDIR / "burst_aware_comparison.png"
    fig.savefig(out, dpi=150)
    print(f"\nfigure -> {out}")
    print(f"total {time.perf_counter()-T0:.0f}s")


if __name__ == "__main__":
    main()
