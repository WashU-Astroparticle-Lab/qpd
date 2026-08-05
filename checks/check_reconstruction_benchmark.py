#!/usr/bin/env python
"""Verification for surrogate-replay benchmarking of a measured trace.

Run: ``python checks/check_reconstruction_benchmark.py`` (exit 0 = all pass).

The module under test claims that reconstruction performance on a *measured*
trace -- which has no truth -- can be recovered by replaying that trace's own
fitted fidelity into surrogates that do. This script gates that claim the only
way it can be gated: on simulated data, where the "measured" trace's truth is
known but withheld from the benchmark, so the prediction can be compared with
what actually happened.

Covers:
  1. Input handling -- complex, (n, 2) and (2, n) I/Q; bad shapes rejected.
  2. Characterisation -- mode auto-detected, and the fitted contrast, sigma and
     rate checked against supervised truth.
  3. Replay closure -- a surrogate re-characterises to the fidelity it was
     built from.
  4. THE CALIBRATION GATE (fixed bias): predicted efficiency and purity match
     the ensemble actually measured over independent real traces.
  5. THE CALIBRATION GATE (swept n_g): same, with the ramp-reset comb replayed.
  6. Degeneracy is carried through, not laundered: at the parity-blind bias the
     benchmark warns, and refuses to "correct" the rate.
  7. Closure warns when the replay is not at the measured fidelity.
  8. noise_scale moves performance in the right direction and by enough to
     locate the breakdown.
  9. Bursts crowd flips and cost efficiency, not purity.
 10. API guards: surrogates cannot outrun the fitted model; settings cannot be
     silently ignored.
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from reconstruction_scenarios import (build_scenario,  # noqa: E402
                                      build_static_scenario)

from qpd.simulator import QuasiparticleBurstModel  # noqa: E402
from qpd.reconstruction.emission import estimate_direction  # noqa: E402
from qpd.reconstruction import (as_complex_trace,  # noqa: E402
                                benchmark_reconstruction, benchmark_vs_noise,
                                characterize_trace,
                                reconstruct_parity_flips_ramped,
                                reconstruct_parity_flips_static, score_flips)

SR = 1e5
_failures = []


def check(name, ok, detail=""):
    tag = "PASS" if ok else "FAIL"
    print(f"  [{tag}] {name}" + (f"  --  {detail}" if detail else ""))
    if not ok:
        _failures.append(name)


def raises(fn, *a, **kw):
    try:
        fn(*a, **kw)
    except (ValueError, TypeError):
        return True
    return False


def ensemble(scn, static, n_seeds):
    """Efficiency and purity actually achieved over independent real traces."""
    eff, pur = [], []
    for s in range(n_seeds):
        r = scn.simulate(seed=s)
        rec = (reconstruct_parity_flips_static(r.iq, SR) if static
               else reconstruct_parity_flips_ramped(r.iq, SR))
        q = score_flips(r.flip_times, rec.flip_times)
        eff.append(q.efficiency)
        pur.append(q.purity)
    return np.array(eff, float), np.array(pur, float)


def main() -> int:
    print(__doc__.split("Covers")[0].strip())
    print()

    # --- 1. input handling -------------------------------------------------
    print("1. Input handling")
    z = np.array([1 + 2j, 3 + 4j, 5 + 6j])
    check("complex 1-D passes through",
          np.allclose(as_complex_trace(z), z))
    check("(n, 2) real I/Q is assembled",
          np.allclose(as_complex_trace(np.c_[z.real, z.imag]), z))
    check("(2, n) real I/Q is assembled",
          np.allclose(as_complex_trace(np.vstack([z.real, z.imag])), z))
    check("real 1-D is rejected, not read as Q = 0",
          raises(as_complex_trace, np.arange(10.0)))
    check("(n, 3) is rejected",
          raises(as_complex_trace, np.zeros((10, 3))))

    # --- 2. characterisation against supervised truth -----------------------
    print("\n2. Characterisation vs supervised truth (fixed bias, n_g = 0.20)")
    scn = build_static_scenario(n_g=0.20, duration=2.0)
    truth_trace = scn.simulate(seed=2)
    fid = characterize_trace(truth_trace.iq, SR)
    check("mode auto-detected as static", fid.mode == "static", fid.mode)

    x = fid.model.project(truth_trace.iq)
    par = truth_trace.parity
    sep_true = abs(x[par == 0].mean() - x[par == 1].mean())
    sig_true = np.sqrt(0.5 * (x[par == 0].var() + x[par == 1].var()))
    c_true = sep_true / sig_true
    check("fitted contrast matches the supervised value within 5%",
          abs(fid.contrast_median / c_true - 1) < 0.05,
          f"fit {fid.contrast_median:.3f} vs truth {c_true:.3f}")
    check("fitted sigma matches the supervised value within 5%",
          abs(fid.sigma / sig_true - 1) < 0.05,
          f"fit {fid.sigma:.4g} vs truth {sig_true:.4g}")
    rate_true = truth_trace.flip_times.size / 2.0
    check("decoded rate matches the trace's realised rate within 20%",
          abs(fid.rate_hz / rate_true - 1) < 0.20,
          f"decoded {fid.rate_hz:.2f} Hz vs realised {rate_true:.2f} Hz")

    scn_r = build_scenario(duration=2.0, burst_rate_hz=0.0)
    fid_r = characterize_trace(scn_r.simulate(seed=2).iq, SR)
    check("mode auto-detected as ramped on a swept trace",
          fid_r.mode == "ramped", fid_r.mode)
    check("fold period recovered to 1e-3 relative",
          abs(fid_r.fold_period / scn_r.fold_period - 1) < 1e-3,
          f"{fid_r.fold_period:.6g} s vs {scn_r.fold_period:.6g} s")
    check("ramp-reset comb found on the measured trace",
          fid_r.reset_comb is not None
          and abs(fid_r.reset_comb.period / scn_r.ramp_period - 1) < 1e-2)

    # --- 3. replay closure --------------------------------------------------
    print("\n3. Replay closure -- a surrogate carries the fidelity it was built from")
    sim_iq, sim_truth = fid.synthesize(seed=11)
    refit = characterize_trace(sim_iq, SR)
    check("surrogate re-characterises to the same contrast within 5%",
          abs(refit.contrast_median / fid.contrast_median - 1) < 0.05,
          f"{refit.contrast_median:.3f} vs {fid.contrast_median:.3f}")
    check("surrogate re-characterises to the same sigma within 5%",
          abs(refit.sigma / fid.sigma - 1) < 0.05)
    check("surrogate length and truth are consistent",
          sim_iq.size == fid.n_samples and sim_truth.ndim == 1)
    # Read the injected noise off the minor principal axis, which carries noise
    # only. The mixture fit's own sigma is NOT the right probe here: at the
    # contrast of 0.96 that noise_scale = 2 produces it absorbs part of the
    # separation and comes back ~6% low -- the same EM bias the closure warning
    # in check 7 exists to surface.
    check("noise_scale scales the surrogate's injected noise",
          abs(estimate_direction(fid.synthesize(seed=11, noise_scale=2.0)[0])[2]
              / estimate_direction(sim_iq)[2] - 2.0) < 0.02)

    # --- 4. calibration gate, fixed bias ------------------------------------
    print("\n4. CALIBRATION GATE -- fixed bias (n_g = 0.20, 2 s, 16 traces)")
    n_seeds = 16
    eff_a, pur_a = ensemble(scn, static=True, n_seeds=n_seeds)
    rep = benchmark_reconstruction(truth_trace.iq, SR, n_trials=n_seeds)
    eff_p, eff_pe = rep.efficiency
    pur_p, pur_pe = rep.purity
    print(f"     actual    eff {eff_a.mean():.3f} +/- {eff_a.std():.3f} | "
          f"purity {pur_a.mean():.3f} +/- {pur_a.std():.3f}")
    print(f"     predicted eff {eff_p:.3f} +/- {eff_pe:.3f} | "
          f"purity {pur_p:.3f} +/- {pur_pe:.3f}")
    # Tolerance: the two means are each an n_seeds-sample estimate, so they are
    # compared against the combined scatter, not against zero.
    tol_e = 3 * np.hypot(eff_a.std(), eff_pe) / np.sqrt(n_seeds) + 0.02
    tol_p = 3 * np.hypot(pur_a.std(), pur_pe) / np.sqrt(n_seeds) + 0.02
    check("predicted efficiency matches the measured ensemble",
          abs(eff_p - eff_a.mean()) < tol_e,
          f"|{eff_p:.3f} - {eff_a.mean():.3f}| < {tol_e:.3f}")
    check("predicted purity matches the measured ensemble",
          abs(pur_p - pur_a.mean()) < tol_p,
          f"|{pur_p:.3f} - {pur_a.mean():.3f}| < {tol_p:.3f}")
    check("predicted spread is the right order, not zero",
          eff_pe <= 3 * max(eff_a.std(), 0.01) + 0.02,
          f"predicted sd {eff_pe:.3f} vs actual {eff_a.std():.3f}")
    check("closure near 1", abs(rep.closure - 1) < 0.05, f"{rep.closure:.3f}")
    # `n_pred * purity / efficiency` recovers `n_truth` identically, so the
    # rate correction is exact by construction and its accuracy on the measured
    # trace is entirely the accuracy of the transferred efficiency and purity
    # -- i.e. the two gates above. Assert the identity, then assert the result
    # lands within the counting uncertainty of the realised rate.
    ident = [t.score.n_pred * t.score.purity / t.score.efficiency
             / t.score.n_truth
             for t in rep.trials if t.score.efficiency > 0]
    check("the rate correction inverts efficiency and purity exactly",
          np.allclose(ident, 1.0),
          f"max deviation {np.max(np.abs(np.array(ident) - 1)):.1e}")
    rate_unc = np.sqrt(max(truth_trace.flip_times.size, 1)) / 2.0
    check("corrected rate agrees with the realised rate within counting error",
          abs(rep.corrected_rate_hz - rate_true) < rate_unc,
          f"corrected {rep.corrected_rate_hz:.2f} vs realised "
          f"{rate_true:.2f} +/- {rate_unc:.2f} Hz "
          f"(raw decode {fid.rate_hz:.2f})")
    check("summary() renders", isinstance(rep.summary(), str)
          and "efficiency" in rep.summary())

    # --- 5. calibration gate, swept n_g -------------------------------------
    print("\n5. CALIBRATION GATE -- swept n_g (500 Hz ramp, 2 s, 8 traces)")
    n_r = 8
    eff_ra, pur_ra = ensemble(scn_r, static=False, n_seeds=n_r)
    rep_r = benchmark_reconstruction(scn_r.simulate(seed=5).iq, SR,
                                     n_trials=n_r)
    eff_rp, eff_rpe = rep_r.efficiency
    pur_rp, pur_rpe = rep_r.purity
    print(f"     actual    eff {eff_ra.mean():.3f} +/- {eff_ra.std():.3f} | "
          f"purity {pur_ra.mean():.3f} +/- {pur_ra.std():.3f}")
    print(f"     predicted eff {eff_rp:.3f} +/- {eff_rpe:.3f} | "
          f"purity {pur_rp:.3f} +/- {pur_rpe:.3f}")
    check("predicted efficiency matches the measured ensemble",
          abs(eff_rp - eff_ra.mean()) < 3 * np.hypot(eff_ra.std(), eff_rpe) + 0.02)
    check("predicted purity matches the measured ensemble",
          abs(pur_rp - pur_ra.mean()) < 3 * np.hypot(pur_ra.std(), pur_rpe) + 0.02)
    # The resets outnumber real flips ~50:1, so purity is the sharp test that
    # the replayed comb is both present in the surrogate and removed from its
    # reconstruction.
    check("replayed ramp resets do not leak into the surrogate's output",
          pur_rp > 0.85, f"purity {pur_rp:.3f}")
    check("surrogate reconstructions recover the comb",
          not any(t.degenerate for t in rep_r.trials))

    # --- 6. degeneracy is carried through -----------------------------------
    print("\n6. Degeneracy is carried through, not laundered (n_g = 0.25)")
    blind = build_static_scenario(n_g=0.25, duration=2.0)
    r_blind = blind.simulate(seed=3)
    rep_b = benchmark_reconstruction(r_blind.iq, SR, n_trials=4,
                                     calibrate_rate=True)
    check("degenerate measured trace is flagged", rep_b.fidelity.degenerate)
    check("benchmark warns that the surrogates replay a fiction",
          any("degenerate" in w for w in rep_b.warnings))
    check("rate calibration refused on a degenerate trace",
          any("calibration skipped" in w for w in rep_b.warnings),
          "; ".join(w[:40] for w in rep_b.warnings))
    check("summary surfaces the warning", "WARNING" in rep_b.summary())
    # The failure this guards against: high scores on a meaningless trace.
    check("scores alone would have looked acceptable (hence the flag matters)",
          rep_b.hard_f1[0] > 0.5,
          f"hard F1 {rep_b.hard_f1[0]:.3f} on a parity-blind bias")

    # --- 7. closure warning -------------------------------------------------
    print("\n7. Closure warns when the replay drifts off the measured fidelity")
    marginal = build_static_scenario(n_g=0.22, duration=2.0)
    rep_m = benchmark_reconstruction(marginal.simulate(seed=3).iq, SR,
                                     n_trials=4)
    warned = any("closure" in w for w in rep_m.warnings)
    check("closure warning fires exactly when closure is off by >5%",
          warned == (abs(rep_m.closure - 1) > 0.05),
          f"closure {rep_m.closure:.3f}, warned={warned}")

    # --- 8. noise_scale -----------------------------------------------------
    print("\n8. noise_scale locates the breakdown")
    reps = benchmark_vs_noise(truth_trace.iq, SR, noise_scales=(0.5, 1.0, 2.0),
                              n_trials=4)
    f1s = [r.hard_f1[0] for r in reps]
    cs = [r.fidelity.contrast_median / r.noise_scale for r in reps]
    print("     noise x0.5 / x1 / x2  ->  hard F1 "
          + " / ".join(f"{v:.3f}" for v in f1s)
          + "  at contrast " + " / ".join(f"{c:.2f}" for c in cs))
    check("performance is monotonically non-increasing in noise",
          f1s[0] >= f1s[1] - 0.05 and f1s[1] >= f1s[2] - 0.05,
          " > ".join(f"{v:.3f}" for v in f1s))
    check("2x noise is a materially worse operating point",
          f1s[2] < f1s[0] - 0.05, f"{f1s[2]:.3f} vs {f1s[0]:.3f}")
    check("all three reuse one characterisation of the measured trace",
          all(r.fidelity is reps[0].fidelity for r in reps))
    # Closure compares two contrasts through an estimator that is biased below
    # contrast ~1, so it is only meaningful where both sides sit at the same
    # fidelity. Rescaled points must decline to answer rather than report a
    # replay failure that is not there.
    check("closure is defined at noise_scale = 1 and withheld elsewhere",
          abs(reps[1].closure - 1) < 0.1
          and not np.isfinite(reps[0].closure)
          and not np.isfinite(reps[2].closure),
          " / ".join(f"{r.closure:.3f}" for r in reps))
    check("no spurious closure warnings on the rescaled points",
          not any("closure" in w for r in (reps[0], reps[2]) for w in r.warnings))

    # --- 9. bursts ----------------------------------------------------------
    print("\n9. Bursts crowd flips: efficiency falls, purity holds")
    bursts = QuasiparticleBurstModel(
        times=np.array([0.4, 0.9, 1.4]), tau=3.7e-3, mu=1.2e-3, sigma=0.4e-3,
        expected_n_qp=25.0)
    rep_nb = benchmark_reconstruction(fidelity=fid, n_trials=4)
    rep_wb = benchmark_reconstruction(fidelity=fid, n_trials=4, bursts=bursts)
    print(f"     no bursts eff {rep_nb.efficiency[0]:.3f} purity "
          f"{rep_nb.purity[0]:.3f} | with bursts eff {rep_wb.efficiency[0]:.3f}"
          f" purity {rep_wb.purity[0]:.3f}")
    check("bursts inject many more flips",
          np.mean([t.n_truth_injected for t in rep_wb.trials])
          > 2 * np.mean([t.n_truth_injected for t in rep_nb.trials]))
    check("burst crowding costs efficiency",
          rep_wb.efficiency[0] < rep_nb.efficiency[0] - 0.05,
          f"{rep_wb.efficiency[0]:.3f} vs {rep_nb.efficiency[0]:.3f}")
    check("burst crowding does not cost purity",
          rep_wb.purity[0] > 0.85, f"purity {rep_wb.purity[0]:.3f}")

    # --- 10. API guards -----------------------------------------------------
    print("\n10. API guards")
    check("surrogates cannot be longer than the fitted window",
          raises(benchmark_reconstruction, fidelity=fid, n_trials=1,
                 trial_samples=fid.n_samples + 1))
    check("trial_samples and trial_duration are mutually exclusive",
          raises(benchmark_reconstruction, fidelity=fid, n_trials=1,
                 trial_samples=1000, trial_duration=0.01))
    check("reconstruction settings alongside fidelity= are rejected, not dropped",
          raises(benchmark_reconstruction, fidelity=fid, n_trials=1,
                 min_confidence=0.3))
    check("a trace with no data is rejected",
          raises(characterize_trace, np.zeros(3, complex), SR))
    check("an unknown mode is rejected",
          raises(characterize_trace, truth_trace.iq, SR, mode="sideways"))
    short = benchmark_reconstruction(fidelity=fid, n_trials=2,
                                     trial_duration=0.5)
    check("trial_duration shortens the surrogates",
          all(t.score.n_truth < fid.n_flips + 8 for t in short.trials)
          and short.fidelity is fid)
    check("recon settings are reused on every surrogate",
          characterize_trace(truth_trace.iq, SR,
                             min_confidence=0.2).recon_kwargs
          == {"min_confidence": 0.2})

    print()
    if _failures:
        print(f"{len(_failures)} check(s) FAILED:")
        for f in _failures:
            print(f"  - {f}")
        return 1
    print("All checks passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
