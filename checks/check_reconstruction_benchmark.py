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
 11. The burst finder: false-burst rate on pure background (and its
     rate-invariance), and real bursts recovered.
 12. The rate and burst-size sweeps, including the multiplicity saturation that
     is the whole point of the burst study.
 13. The three diagnostic figures render.
 14b. Error-bar semantics: scatter vs precision-of-estimate, and that pooling
      weights trials by evidence rather than one-vote-each.
 14. The two decoding rules: identical at good contrast, and diverging in the
     documented direction (recall up, purity down) near the noise floor.
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from reconstruction_scenarios import (build_scenario,  # noqa: E402
                                      build_static_scenario)

from qpd.simulator import QuasiparticleBurstModel  # noqa: E402
from qpd.simulator.parity import generate_parity_trajectory  # noqa: E402
from qpd.reconstruction.emission import estimate_direction  # noqa: E402
from qpd.reconstruction import (as_complex_trace,  # noqa: E402
                                benchmark_reconstruction, benchmark_vs_noise,
                                characterize_trace, detect_bursts,
                                match_bursts, plot_burst_efficiency,
                                plot_burst_multiplicity,
                                plot_efficiency_vs_rate,
                                reconstruct_parity_flips_ramped,
                                reconstruct_parity_flips_static, score_flips,
                                sweep_burst_size, sweep_rate)

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

    # Regression: a fixed bias must not be auto-detected as swept. The
    # telegraph's own Lorentzian spectrum peaks at the bottom of the period
    # search band, and on a longer trace it is reported as a highly significant
    # period of order the trace length. Following it sends the trace down the
    # ramped pipeline, which models telegraph noise as a ramp and inflates the
    # decoded rate several-fold -- silently. 5 s is where this first bit.
    for ng in (0.0, 0.10, 0.20):
        for dur in (2.0, 5.0):
            f = characterize_trace(
                build_static_scenario(n_g=ng, duration=dur).simulate(seed=1).iq,
                SR)
            check(f"fixed bias n_g={ng} over {dur:.0f} s stays static",
                  f.mode == "static", f"mode={f.mode}, rate={f.rate_hz:.1f} Hz")

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

    # --- 11. burst finder ---------------------------------------------------
    print("\n11. Burst finder on a flip train")
    # The number that matters: how often pure background fakes a burst. It must
    # also be rate-invariant, since the linking distance scales with the rate.
    for rate in (5.0, 50.0, 500.0):
        n_false, n_tr, D = 0, 120, 5.0
        for s in range(n_tr):
            r = np.random.default_rng(s)
            t = np.sort(r.uniform(0.0, D, r.poisson(rate * D)))
            n_false += len(detect_bursts(t, rate, duration=D))
        check(f"false bursts on pure {rate:.0f} Hz background stay rare",
              n_false / n_tr < 0.05, f"{n_false / n_tr:.3f} per {D:.0f} s trace")

    # And that it finds real ones. On a truth train (no reconstruction), a
    # 20-quasiparticle burst is unmistakable.
    D, rate = 5.0, 10.0
    grid = np.arange(int(D * SR)) / SR
    found, total = 0, 0
    for s in range(8):
        r = np.random.default_rng(4000 + s)
        bm = QuasiparticleBurstModel(times=np.arange(0.5, D - 0.5, 0.5),
                                     tau=3.7e-3, mu=1.2e-3, sigma=0.4e-3,
                                     expected_n_qp=20.0)
        ev, bt = bm.sample(r)
        _, flips = generate_parity_trajectory(grid, rate, rate, r,
                                              extra_flip_times=ev,
                                              return_flip_times=True)
        for m in match_bursts(bt, detect_bursts(flips, rate, duration=D)):
            total += 1
            found += m.detected
    check("20-qp bursts are found on a truth train",
          found / total > 0.95, f"{found}/{total}")
    check("match_bursts drops empty (n_qp = 0) truth bursts",
          all(m.n_qp_true > 0 for m in match_bursts(
              [type("B", (), {"n_qp": 0, "t_start": 0.0, "t_end": 0.1})()],
              [])) is True)
    check("a zero background rate is rejected, not silently allowed",
          raises(detect_bursts, np.array([0.0, 0.001, 0.002]), 0.0))
    check("min_flips below 2 is rejected",
          raises(detect_bursts, np.array([0.0, 0.001]), 10.0, min_flips=1))

    # --- 12. sweeps ---------------------------------------------------------
    print("\n12. Rate and burst-size sweeps")
    fid5 = characterize_trace(
        build_static_scenario(n_g=0.0, duration=5.0).simulate(seed=1).iq, SR)
    reports = sweep_rate(fid5, [10.0, 1e3, 1e4], n_trials=3)
    effs = [r.efficiency[0] for r in reports]
    print("     rate 10 / 1k / 10k Hz -> efficiency "
          + " / ".join(f"{e:.3f}" for e in effs))
    check("flip efficiency does not increase with background rate",
          effs[0] >= effs[1] - 0.02 and effs[1] >= effs[2] - 0.02,
          " > ".join(f"{e:.3f}" for e in effs))
    check("crowding costs recall before precision",
          reports[-1].purity[0] > reports[-1].efficiency[0] - 1e-9,
          f"purity {reports[-1].purity[0]:.3f} vs "
          f"efficiency {reports[-1].efficiency[0]:.3f}")
    check("the rate is pinned, not jittered, across a rate sweep",
          all(not r.rate_jitter for r in reports))
    # The tolerance varies along the sweep and the figure deliberately does not
    # say so, so it has to stay recoverable from the reports themselves.
    tols = [r.tol for r in reports]
    check("each point records the tolerance it was scored at",
          tols[0] > tols[1] > tols[2] and all(t > 0 for t in tols),
          " > ".join(f"{t * 1e6:.1f}us" for t in tols))
    check("adaptive_tol=False scores every rate at one fixed window",
          len({r.tol for r in sweep_rate(fid5, [10.0, 100.0], n_trials=1,
                                         adaptive_tol=False)}) == 1)

    pts = sweep_burst_size(fid5, [3, 8, 30, 80], n_trials=3)
    print("     n_qp 3 / 8 / 30 / 80 -> burst efficiency "
          + " / ".join(f"{p.efficiency:.2f}" for p in pts)
          + "  bias " + " / ".join(f"{p.bias:+.1f}" for p in pts))
    check("burst detection efficiency rises with multiplicity",
          pts[0].efficiency < pts[1].efficiency <= pts[2].efficiency,
          " < ".join(f"{p.efficiency:.3f}" for p in pts[:3]))
    check("large bursts are always found",
          pts[-1].efficiency > 0.95, f"{pts[-1].efficiency:.3f}")
    # The physics this whole study exists to expose: unresolved pairs cancel in
    # the parity, so recovered multiplicity falls progressively short.
    check("multiplicity saturates -- big bursts are under-counted",
          pts[-1].bias < -3.0, f"bias {pts[-1].bias:+.1f} at n_qp = 80")
    check("small bursts are not under-counted the same way",
          pts[0].bias > pts[-1].bias, f"{pts[0].bias:+.1f} vs {pts[-1].bias:+.1f}")
    check("the burst study runs at the background rate it was given",
          all(abs(p.background_rate_hz - fid5.rate_hz) < 1e-9 for p in pts))

    # --- 13. diagnostic plots ----------------------------------------------
    print("\n13. Diagnostic plots render")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    try:
        f1, _ = plot_efficiency_vs_rate(reports)
        f2, _ = plot_burst_efficiency(pts)
        f3, _ = plot_burst_multiplicity(pts)
        check("all three diagnostic figures render",
              all(f is not None for f in (f1, f2, f3)))
        plt.close("all")
    except Exception as exc:  # noqa: BLE001
        check("all three diagnostic figures render", False, repr(exc))
    check("an empty sweep is rejected rather than plotted",
          raises(plot_efficiency_vs_rate, []))

    # --- 14b. error-bar semantics -------------------------------------------
    print("\n14b. Error bars: scatter vs precision of the estimate")
    fe = characterize_trace(
        build_static_scenario(n_g=0.18, duration=5.0, sample_rate=1e4,
                              tunnel_rate_hz=66.0).simulate(seed=1).iq, 1e4)
    small = sweep_rate(fe, [3.0], n_trials=30, seed=3)[0]
    big = sweep_rate(fe, [3.0], n_trials=120, seed=3)[0]
    # The trial-to-trial SD estimates a population scatter, so 4x the trials
    # must NOT shrink it appreciably; the pooled/bootstrap errors must.
    check("trial scatter (sd) does not shrink with n_trials",
          big.purity[1] > 0.5 * small.purity[1],
          f"{small.purity[1]:.4f} -> {big.purity[1]:.4f}")
    check("pooled binomial error shrinks with n_trials",
          big.purity_pooled[1] < small.purity_pooled[1],
          f"{small.purity_pooled[1]:.4f} -> {big.purity_pooled[1]:.4f}")
    check("bootstrap error shrinks with n_trials",
          big.purity_bootstrap[1] < small.purity_bootstrap[1],
          f"{small.purity_bootstrap[1]:.4f} -> {big.purity_bootstrap[1]:.4f}")
    check("bootstrap is never tighter than binomial (it adds clustering)",
          big.purity_bootstrap[1] >= 0.9 * big.purity_pooled[1],
          f"boot {big.purity_bootstrap[1]:.4f} vs binom {big.purity_pooled[1]:.4f}")
    check("pooled and bootstrap agree on the point estimate",
          abs(big.purity_pooled[0] - big.purity_bootstrap[0]) < 1e-12)
    check("pooled F1 is consistent with pooled precision and recall",
          abs(big.hard_f1_pooled
              - 2 * big.efficiency_pooled[0] * big.purity_pooled[0]
              / (big.efficiency_pooled[0] + big.purity_pooled[0])) < 1e-12)

    # Pooling must weight by evidence: a trial with no truth and a burst of
    # spurious flips has to count for more than one vote among many.
    n_pred = np.array([t.score.n_pred for t in big.trials])
    if n_pred.max() > 3 * max(np.median(n_pred), 1):
        check("pooling weights trials by evidence, mean-of-ratios does not",
              big.purity_pooled[0] <= big.purity[0] + 1e-9,
              f"pooled {big.purity_pooled[0]:.3f} vs mean-of-ratios "
              f"{big.purity[0]:.3f} (max n_pred {n_pred.max()}, "
              f"median {int(np.median(n_pred))})")
    check("an unknown err= is rejected",
          raises(plot_efficiency_vs_rate, [big], err="nope"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as _plt
    for kind in ("bootstrap", "binomial", "sd"):
        plot_efficiency_vs_rate([big], err=kind)
    _plt.close("all")
    check("all three err= modes render", True)

    # --- 14. decoding rule ---------------------------------------------------
    print("\n14. Viterbi vs forward-backward decoding")
    check("an unknown decoder is rejected",
          raises(reconstruct_parity_flips_static, truth_trace.iq, SR,
                 decoder="nope"))
    # Agreement is only expected where the evidence is strong, so this uses the
    # best-contrast bias (n_g = 0). At the n_g = 0.20 trace above the contrast
    # is 1.18 -- already inside the band where the two rules legitimately
    # differ, which is the point of the divergence test further down.
    good = build_static_scenario(n_g=0.0, duration=5.0).simulate(seed=1)
    rv = reconstruct_parity_flips_static(good.iq, SR, decoder="viterbi")
    rp = reconstruct_parity_flips_static(good.iq, SR, decoder="posterior")
    check("both decoders run and record which was used",
          rv.diagnostics["decoder"] == "viterbi"
          and rp.diagnostics["decoder"] == "posterior")
    # The rate is estimated from the posterior either way, so switching the
    # decoder must not move it -- that is what makes the comparison a clean
    # test of the decoding step alone.
    check("the decoder does not change the fitted rate or contrast",
          rv.rate_hz == rp.rate_hz and rv.contrast == rp.contrast,
          f"{rv.rate_hz:.4f} vs {rp.rate_hz:.4f} Hz")
    check("the posterior path is the marginal thresholded at 1/2",
          np.array_equal(rp.branch, (rp.posterior > 0.5).astype(np.int8)))
    check("at good contrast the two rules return the same sequence",
          np.array_equal(rv.branch, rp.branch),
          f"contrast {rv.contrast:.2f}")

    # They must diverge in the documented direction once the evidence is
    # marginal: the posterior rule follows short excursions, buying recall and
    # paying purity. This is the claim §12e and the notebook rest on.
    fnoisy = characterize_trace(
        build_static_scenario(n_g=0.0, duration=5.0).simulate(seed=1).iq, SR)
    lo = {d: benchmark_reconstruction(
              fidelity=characterize_trace(
                  build_static_scenario(n_g=0.0, duration=5.0).simulate(seed=1).iq,
                  SR, decoder=d),
              n_trials=6, noise_scale=5.6, rate_jitter=False)
          for d in ("viterbi", "posterior")}
    ev, pv = lo["viterbi"].efficiency[0], lo["viterbi"].purity[0]
    ep, pp = lo["posterior"].efficiency[0], lo["posterior"].purity[0]
    print(f"     at contrast {fnoisy.contrast_median / 5.6:.2f}: "
          f"viterbi {ev:.3f}/{pv:.3f}, posterior {ep:.3f}/{pp:.3f} (eff/pur)")
    check("near the noise floor the posterior rule has the higher recall",
          ep >= ev - 1e-9, f"{ep:.3f} vs {ev:.3f}")
    check("...and pays for it in purity",
          pp < pv, f"{pp:.3f} vs {pv:.3f}")
    check("...so Viterbi wins on F1 there",
          lo["viterbi"].hard_f1[0] > lo["posterior"].hard_f1[0],
          f"{lo['viterbi'].hard_f1[0]:.3f} vs {lo['posterior'].hard_f1[0]:.3f}")
    check("the decoder reaches every surrogate through recon_kwargs",
          all(t.score.n_truth >= 0 for t in lo["posterior"].trials)
          and lo["posterior"].fidelity.recon_kwargs == {"decoder": "posterior"})

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
