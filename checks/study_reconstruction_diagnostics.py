#!/usr/bin/env python
"""Diagnostic plots for a measured trace: rate response and burst response.

Run: ``python checks/study_reconstruction_diagnostics.py [rate|burst|all]``

Everything is anchored on one "measured" trace -- here a simulated one, but the
entry point is the same for real data. The trace is characterised once, and its
fitted fidelity drives every sweep, so all three figures describe the readout
you actually have rather than a generic device.

Figures written to ``docs/figures/``:

``efficiency_vs_rate.png``
    Flip detection efficiency and purity against the Poisson background
    tunnelling rate.

``burst_efficiency.png``
    Burst-level detection efficiency against burst quasiparticle number, at the
    background rate measured from the trace.

``burst_multiplicity.png``
    True against reconstructed burst multiplicity, and the bias.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib  # noqa: E402

matplotlib.use("Agg")

from reconstruction_scenarios import build_static_scenario  # noqa: E402

from qpd.reconstruction import (benchmark_reconstruction,  # noqa: E402
                                characterize_trace, plot_burst_efficiency,
                                plot_burst_multiplicity,
                                plot_efficiency_vs_rate, sweep_burst_size,
                                sweep_rate)

FIGDIR = Path(__file__).resolve().parents[1] / "docs" / "figures"
SR = 1e5
RATES = [1, 3, 10, 30, 100, 300, 1e3, 3e3, 1e4, 3e4]
N_QP = [2, 3, 4, 5, 6, 8, 10, 12, 16, 20, 30, 50, 80, 120]


def measured_trace():
    """The stand-in for the user's data: one 5 s fixed-bias trace."""
    scn = build_static_scenario(n_g=0.0, duration=5.0, tunnel_rate_hz=10.0)
    return scn.simulate(seed=1)


def main(which="all") -> int:
    FIGDIR.mkdir(parents=True, exist_ok=True)
    result = measured_trace()

    fid = characterize_trace(result.iq, SR)
    print(f"Measured trace: mode={fid.mode}, contrast={fid.contrast_median:.2f}, "
          f"decoded rate={fid.rate_hz:.2f} Hz, {fid.n_flips} flips")

    # The background rate the burst study is run at comes from the data, and
    # from the de-biased value rather than the raw decode.
    rep = benchmark_reconstruction(fidelity=fid, n_trials=8)
    bg = rep.corrected_rate_hz
    print(f"  efficiency {rep.efficiency[0]:.3f} +/- {rep.efficiency[1]:.3f}, "
          f"corrected background rate {bg:.2f} Hz "
          f"(true 10.0 Hz, realised {result.flip_times.size / 5:.2f} Hz)")

    if which in ("rate", "all"):
        print("\n-- flip efficiency vs background rate --")
        reports = sweep_rate(fid, RATES, n_trials=5)
        print(f"  {'rate/Hz':>9s} {'effic':>7s} {'purity':>7s} {'F1':>7s}")
        for r, rt in zip(reports, RATES):
            print(f"  {rt:9.4g} {r.efficiency[0]:7.3f} "
                  f"{r.purity[0]:7.3f} {r.hard_f1[0]:7.3f}")
        fig, _ = plot_efficiency_vs_rate(reports)
        fig.savefig(FIGDIR / "efficiency_vs_rate.png", dpi=200)
        print(f"  wrote {FIGDIR / 'efficiency_vs_rate.png'}")

    if which in ("burst", "all"):
        print("\n-- burst detection and multiplicity vs burst size --")
        pts = sweep_burst_size(fid, N_QP, background_rate_hz=bg, n_trials=6)
        print(f"  {'n_qp':>5s} {'bursts':>7s} {'burst_eff':>18s} "
              f"{'<true>':>7s} {'<det>':>7s} {'bias':>8s}")
        for p in pts:
            print(f"  {p.n_qp_expected:5.0f} {p.n_bursts:7d} "
                  f"{p.efficiency:9.3f} +/- {p.efficiency_err:.3f} "
                  f"{p.mean_n_qp_true:7.1f} {p.mean_n_qp_detected:7.1f} "
                  f"{p.bias:+8.1f}")
        fig, _ = plot_burst_efficiency(pts)
        fig.savefig(FIGDIR / "burst_efficiency.png", dpi=200)
        print(f"  wrote {FIGDIR / 'burst_efficiency.png'}")
        fig, _ = plot_burst_multiplicity(pts)
        fig.savefig(FIGDIR / "burst_multiplicity.png", dpi=200)
        print(f"  wrote {FIGDIR / 'burst_multiplicity.png'}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1] if len(sys.argv) > 1 else "all"))
