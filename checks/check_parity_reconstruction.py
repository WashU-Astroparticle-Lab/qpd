#!/usr/bin/env python
"""Verification for blind parity-flip timing reconstruction.

Run: ``python checks/check_parity_reconstruction.py`` (exit 0 = all pass).

Covers, on the 500 Hz-ramp reference scenario:
  1. Nuisance calibration -- fold period and ramp-reset comb recovered blind.
  2. Clean telegraph (10 Hz tunnelling): near-perfect detection and timing.
  3. Quasiparticle bursts: efficiency drops (unresolvable close pairs) but
     purity and timing stay intact.
  4. Offset-charge jumps as nuisances, including the worst case delta = 0.25,
     which blends the splitting profile with its own half-period shift.
  5. Ramp resets never leak into the output: with 50x more resets than real
     flips, purity is the sharp test.
  6. A second seed, so nothing above is seed-tuned.
  7. Plotting helpers render, and reject an out-of-range window.
  8. The constant-bias entry point: two stationary blobs, and the parity-blind
     bias point flagged rather than silently returning nonsense.
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from reconstruction_scenarios import (build_scenario,  # noqa: E402
                                      build_static_scenario)

from qpd.reconstruction import (reconstruct_parity_flips_ramped,  # noqa: E402
                               reconstruct_parity_flips_static, score_flips)

_failures = []


def check(name, ok, detail=""):
    tag = "PASS" if ok else "FAIL"
    print(f"  [{tag}] {name}" + (f"  --  {detail}" if detail else ""))
    if not ok:
        _failures.append(name)


def run(scn, seed=0):
    result = scn.simulate(seed=seed)
    rec = reconstruct_parity_flips_ramped(result.iq, scn.sim.sample_rate)
    return result, rec, score_flips(result.flip_times, rec.flip_times)


def main() -> int:
    print(__doc__.split("Covers")[0].strip())
    print()

    # --- 1. blind nuisance calibration -------------------------------------
    print("1. Blind nuisance calibration (500 Hz ramp, 10 Hz tunnelling)")
    scn = build_scenario(duration=5.0, burst_rate_hz=0.0)
    res, rec, score = run(scn)
    p_err = abs(rec.emission.fold_period / scn.fold_period - 1.0)
    check("fold period recovered to <1e-4 relative", p_err < 1e-4,
          f"rel err {p_err:.2e} (P = {scn.fold_period * 1e6:.3f} us)")
    comb = rec.reset_comb
    check("ramp-reset comb detected", comb is not None)
    if comb is not None:
        t_err = abs(comb.period / scn.ramp_period - 1.0)
        check("reset comb spans the right number of half-pairs",
              comb.n_fold == int(scn.half_pairs_per_ramp),
              f"n_fold = {comb.n_fold} (expected {int(scn.half_pairs_per_ramp)})")
        check("reset period recovered to <1e-3 relative", t_err < 1e-3,
              f"rel err {t_err:.2e} (T = {comb.period * 1e3:.4f} ms)")
    rate_err = abs(rec.rate_hz / 10.0 - 1.0)
    check("tunnelling rate recovered within 40%", rate_err < 0.4,
          f"{rec.rate_hz:.1f} Hz vs 10 Hz truth")

    # --- 2. clean telegraph -------------------------------------------------
    print("\n2. Clean telegraph, no bursts, no jumps")
    check("hard F1 >= 0.95", score.hard_f1 >= 0.95, f"F1 = {score.hard_f1:.3f}")
    check("purity >= 0.97", score.purity >= 0.97, f"purity = {score.purity:.3f}")
    check("timing rms <= 50 us", score.rms_s <= 50e-6,
          f"rms = {score.rms_s * 1e6:.1f} us, bias = {score.bias_s * 1e6:+.1f} us")
    n_reset = int(round(scn.duration / scn.ramp_period))
    check("false positives far below the reset count",
          score.n_pred - score.n_matched < 0.01 * n_reset,
          f"{score.n_pred - score.n_matched} false vs {n_reset} resets in trace")

    # --- 3. quasiparticle bursts -------------------------------------------
    print("\n3. With quasiparticle bursts (1 Hz onsets, 15 tunnels each)")
    scn_b = build_scenario(duration=5.0, burst_rate_hz=1.0, burst_n_qp=15.0)
    _, rec_b, score_b = run(scn_b)
    check("hard F1 >= 0.88", score_b.hard_f1 >= 0.88,
          f"F1 = {score_b.hard_f1:.3f}, efficiency = {score_b.efficiency:.3f}")
    check("purity >= 0.95 (bursts must not manufacture events)",
          score_b.purity >= 0.95, f"purity = {score_b.purity:.3f}")
    check("timing rms <= 50 us", score_b.rms_s <= 50e-6,
          f"rms = {score_b.rms_s * 1e6:.1f} us")

    # --- 4. offset-charge jumps as nuisances -------------------------------
    print("\n4. Offset-charge jumps (nuisances; amplitudes not reconstructed)")
    for delta in (0.125, 0.25, 0.40):
        scn_j = build_scenario(duration=5.0, burst_rate_hz=1.0,
                               charge_jumps=(np.array([2.5]),
                                             np.array([delta])))
        _, rec_j, score_j = run(scn_j)
        found = rec_j.charge_jump_times
        near = (np.min(np.abs(found - 2.5)) if found.size else np.inf)
        check(f"delta = {delta}: hard F1 >= 0.88", score_j.hard_f1 >= 0.88,
              f"F1 = {score_j.hard_f1:.3f}, purity = {score_j.purity:.3f}")
        check(f"delta = {delta}: jump located within 1 ms", near < 1e-3,
              f"nearest detected jump {near * 1e3:.3f} ms away "
              f"({found.size} found)")

    # --- 5. second seed -----------------------------------------------------
    print("\n5. Second seed (guards against seed-tuning)")
    scn_s = build_scenario(duration=5.0, burst_rate_hz=1.0, seed=7)
    _, rec_s, score_s = run(scn_s, seed=7)
    check("hard F1 >= 0.88 on seed 7", score_s.hard_f1 >= 0.88,
          f"F1 = {score_s.hard_f1:.3f}, purity = {score_s.purity:.3f}")
    check("reset comb detected on seed 7", rec_s.reset_comb is not None)

    # --- 6. constant-bias reconstruction -----------------------------------
    print("\n6. Constant bias (reconstruct_parity_flips_static)")
    scn_c = build_static_scenario(n_g=0.0, duration=5.0, tunnel_rate_hz=10.0)
    res_c = scn_c.simulate(seed=0)
    rec_c = reconstruct_parity_flips_static(res_c.iq, scn_c.sim.sample_rate)
    score_c = score_flips(res_c.flip_times, rec_c.flip_times)
    check("n_g=0: hard F1 >= 0.95", score_c.hard_f1 >= 0.95,
          f"F1 = {score_c.hard_f1:.3f}, contrast = {rec_c.contrast:.2f}")
    check("n_g=0: timing rms <= 20 us", score_c.rms_s <= 20e-6,
          f"rms = {score_c.rms_s * 1e6:.1f} us")
    check("n_g=0: rate recovered within 40%", abs(rec_c.rate_hz / 10.0 - 1) < 0.4,
          f"{rec_c.rate_hz:.1f} Hz vs 10 Hz truth")
    check("n_g=0: not flagged degenerate", not rec_c.degenerate,
          f"detectability = {rec_c.diagnostics['detectability']:.0f}")

    # The parity-blind charge: the two branches coincide, so this is
    # unrecoverable in principle. What matters is that it is *flagged* -- the
    # fitted contrast alone would not reveal it, since EM splits a single blob
    # into a spurious pair with contrast ~0.9 either way.
    scn_d = build_static_scenario(n_g=0.25, duration=5.0, tunnel_rate_hz=10.0)
    res_d = scn_d.simulate(seed=0)
    rec_d = reconstruct_parity_flips_static(res_d.iq, scn_d.sim.sample_rate)
    check("n_g=0.25 (parity-blind) flagged degenerate", rec_d.degenerate,
          f"detectability = {rec_d.diagnostics['detectability']:.1f}, "
          f"fitted contrast = {rec_d.contrast:.2f} (spurious)")

    # Cross-check: the static routine must not be used on a ramped trace.
    scn_r = build_scenario(duration=3.0, burst_rate_hz=0.0)
    res_r = scn_r.simulate(seed=0)
    rec_r = reconstruct_parity_flips_static(res_r.iq, scn_r.sim.sample_rate)
    score_r = score_flips(res_r.flip_times, rec_r.flip_times)
    check("static routine on a ramped trace does not silently look fine",
          rec_r.degenerate or score_r.hard_f1 < 0.5,
          f"F1 = {score_r.hard_f1:.3f}, degenerate = {rec_r.degenerate}")

    # --- 7. plotting helpers -----------------------------------------------
    print("\n7. Plotting helpers (headless smoke test)")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from qpd.reconstruction import plot_iq_plane, plot_trace_with_flips

    try:
        fig, axs = plot_trace_with_flips(
            res_c.iq, scn_c.sim.sample_rate, rec_c.flip_times,
            truth_times=res_c.flip_times, window=(0.0, 0.3), smooth_hz=300.0)
        ok_raw = len(np.atleast_1d(axs)) == 2
        plt.close(fig)
        fig, ax = plot_trace_with_flips(
            res.iq, scn.sim.sample_rate, rec.flip_times,
            truth_times=res.flip_times, window=(0.010, 0.013),
            projected=True, emission=rec.emission,
            confidence=rec.confidence)
        plt.close(fig)
        fig, ax = plot_iq_plane(res_c.iq, branch=rec_c.branch, model=rec_c.model)
        plt.close(fig)
        check("plot_trace_with_flips / plot_iq_plane render", ok_raw,
              "raw I/Q (2 panels), projected, and I/Q plane all rendered")
    except Exception as exc:  # noqa: BLE001 - a render failure is the finding
        check("plot_trace_with_flips / plot_iq_plane render", False, repr(exc))

    # A window outside the trace is a caller error, not a silent empty plot.
    try:
        plot_trace_with_flips(res_c.iq, scn_c.sim.sample_rate,
                              window=(99.0, 99.1))
        check("out-of-range window raises", False, "no exception raised")
    except ValueError:
        check("out-of-range window raises", True)
    except Exception as exc:  # noqa: BLE001
        check("out-of-range window raises", False, f"wrong type: {exc!r}")

    print()
    if _failures:
        print(f"RESULT: {len(_failures)} check(s) FAILED: {_failures}")
        return 1
    print("RESULT: all checks PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
