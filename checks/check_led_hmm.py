#!/usr/bin/env python
"""Verify the LED-window HMM decoder and the blind LED-phase estimator.

Covers ``qpd.mlebench.led_hmm`` (``predict_led_hmm``, ``predict_led_hmm_blind``)
and ``qpd.mlebench.led_phase`` (``estimate_led_phase``) on records from the
LED-pulsing simulator, with an inexpensive charge-dependent dispersive shift.

Run: ``python checks/check_led_hmm.py`` (exit 0 = all checks pass). Seeds are
fixed, and the thresholds sit well outside the spread between seeds.
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from qpd import QPD  # noqa: E402
from qpd.mlebench.grade import match_times  # noqa: E402
from qpd.mlebench.led_hmm import predict_led_hmm, predict_led_hmm_blind  # noqa: E402
from qpd.mlebench.led_phase import estimate_led_phase  # noqa: E402
from qpd.reconstruction.static_bias import reconstruct_parity_flips_static  # noqa: E402
from qpd.simulator import (  # noqa: E402
    LEDPulseModel,
    ResonatorConfig,
    VNASimulator,
    WhiteGaussianNoise,
    notch_s21,
)

_failures = []
FS, PERIOD, TOL = 100_000.0, 0.30005, 200e-6
SF = dict(pin_quiet_rate=True, smooth_fit_s=2e-3, min_smoothed_contrast=2.5)
RES = ResonatorConfig(f_r=6.936e9, q_i=1e5, q_c_abs=1.2e5)
SHIFT_HZ = 3e4


def check(name, ok, detail=""):
    tag = "PASS" if ok else "FAIL"
    print(f"  [{tag}] {name}" + (f"  --  {detail}" if detail else ""))
    if not ok:
        _failures.append(name)


def rejects(name, action):
    try:
        action()
    except ValueError:
        check(name, True)
    else:
        check(name, False)


class FastSimulator(VNASimulator):
    """Parity shifts of +-SHIFT_HZ * cos(2 pi n_g): parity-blind charges at n_g = +-1/4."""

    def _chi_grid(self, n_g_lo, n_g_hi):
        grid = np.linspace(n_g_lo - 0.01, n_g_hi + 0.01, 401)
        shift = SHIFT_HZ * np.cos(2 * np.pi * grid)
        return grid, shift, -shift


def max_separation():
    shift = SHIFT_HZ * np.cos(2 * np.pi * np.linspace(-0.5, 0.5, 2001))
    return float(np.max(np.abs(notch_s21(RES.f_r, RES.f_r + shift, RES)
                               - notch_s21(RES.f_r, RES.f_r - shift, RES))))


def simulate(seed, snr_max, duration=2.0, n_qp=3.0):
    """One record: (iq, true flip times, true first LED onset). snr_max = best-charge separation / noise sigma."""
    first = float(np.random.default_rng(seed).uniform(0, PERIOD))
    sim = FastSimulator(
        qpd=QPD(e_j_hz=6.95e9, e_c_hz=0.695e9, coupling_g_hz=40e6), resonator=RES,
        f_drive=RES.f_r, sample_rate=FS, noise=WhiteGaussianNoise(sigma=max_separation() / snr_max),
        led_pulsing=LEDPulseModel(first_pulse_time=first, expected_n_qp=n_qp, charge_jump_rate_hz=0.0))
    res = sim.simulate(duration=duration, seed=seed)
    return res.iq, res.flip_times, first


def veto(times, first, window_s=300e-6):
    """Drop detections in [onset, onset + window_s) after each onset of the schedule starting at ``first``."""
    times = np.asarray(times, float)
    since = (times - first) % PERIOD
    return times[~((times >= first - 1e-12) & (since < window_s))] if times.size else times


def counts(truth, pred):
    m = match_times(truth, pred, TOL)
    return m.tp, len(pred) - m.tp, len(truth) - m.tp


def f1(tp, fp, fn):
    return 2 * tp / max(2 * tp + fp + fn, 1)


def summary(tp, fp, fn):
    """Efficiency (matched / true), purity (matched / reported) and F1 of pooled counts."""
    return (f"eff {tp / max(tp + fn, 1):.3f}, pur {tp / max(tp + fp, 1):.3f}, F1 {f1(tp, fp, fn):.3f}, "
            f"FP {fp}, FN {fn}")


def phase_error(est, true):
    return ((est - true + PERIOD / 2) % PERIOD) - PERIOD / 2


def main():
    print("== contract")
    iq, truth, first = simulate(11, 3.0, duration=1.0)
    rejects("dead_detectability <= 0 rejected",
            lambda: predict_led_hmm(iq, FS, first_pulse_time_s=first, dead_detectability=0))
    rejects("smooth_fit_s without pin_quiet_rate rejected",
            lambda: predict_led_hmm(iq, FS, first_pulse_time_s=first, smooth_fit_s=2e-3))
    rejects("min_smoothed_contrast without smooth_fit_s rejected",
            lambda: predict_led_hmm(iq, FS, first_pulse_time_s=first, pin_quiet_rate=True, min_smoothed_contrast=2.5))
    rejects("unknown phase method rejected", lambda: estimate_led_phase(iq, FS, PERIOD, method="bogus"))

    print("== default decoder = independent per-window static decode")
    ev, info = predict_led_hmm(iq, FS, first_pulse_time_s=first)
    edges = np.unique(np.r_[0, np.searchsorted(np.arange(iq.size) / FS, np.arange(first, iq.size / FS, PERIOD)), iq.size])
    manual = np.concatenate([reconstruct_parity_flips_static(
        iq[a:b], FS, t0=a / FS, p_flip_init=20.0 / FS, n_rate_iterations=4, min_confidence=0.0, segment_blocks=1,
        decoder="viterbi", burst_aware=True, p_burst=0.3, burst_rate_hz=1.0 / PERIOD, burst_tau=3.7e-3).flip_times
        for a, b in zip(edges[:-1], edges[1:]) if b - a >= 8])
    check("default output = manual per-window decode", np.array_equal(np.sort(ev), np.sort(manual)),
          f"{len(ev)} flips, {info['segments']} windows")

    print("== HMM-SF, known phase")
    for snr, seeds, min_f1, fp_ratio in ((3.0, range(20, 26), 0.80, 0.5), (1.0, range(30, 36), None, 0.2)):
        tot = {"plain": np.zeros(3, int), "sf": np.zeros(3, int)}
        for s in seeds:
            iq, truth, first = simulate(s, snr)
            tot["plain"] += counts(truth, predict_led_hmm(iq, FS, first_pulse_time_s=first)[0])
            tot["sf"] += counts(truth, predict_led_hmm(iq, FS, first_pulse_time_s=first, **SF)[0])
        p, q = tot["plain"], tot["sf"]
        print(f"    SNR_max {snr:g}  plain:  {summary(*p)}\n    SNR_max {snr:g}  HMM-SF: {summary(*q)}")
        detail = f"F1 plain {f1(*p):.3f} -> HMM-SF {f1(*q):.3f}; false positives {p[1]} -> {q[1]}"
        if min_f1 is not None:
            check(f"SNR_max {snr:g}: HMM-SF F1 >= {min_f1}", f1(*q) >= min_f1, detail)
        check(f"SNR_max {snr:g}: HMM-SF false positives <= {fp_ratio} x plain", q[1] <= fp_ratio * p[1], detail)

    print("== blind LED phase")
    errs, f1s = [], {"known": np.zeros(3, int), "blind": np.zeros(3, int)}
    for s in range(40, 46):
        iq, truth, first = simulate(s, 3.0)
        est, _ = estimate_led_phase(iq, FS, PERIOD)
        errs.append(phase_error(est, first))
        v1, _ = estimate_led_phase(iq, FS, PERIOD, method="v1")
        check(f"seed {s}: 'v1' estimate lies in [0, period)", 0 <= v1 < PERIOD, f"{v1:.6f} s")
        ev, info = predict_led_hmm_blind(iq, FS, **SF)
        check(f"seed {s}: blind schedule = estimate - 100 us",
              abs(phase_error(info["first_pulse_time_s"], info["estimated_onset_s"] - 100e-6)) < 1e-12)
        # Each decoder's 300 us onset veto sits at the schedule it used.
        f1s["blind"] += counts(truth, veto(ev, info["first_pulse_time_s"]))
        f1s["known"] += counts(truth, veto(predict_led_hmm(iq, FS, first_pulse_time_s=first, **SF)[0], first))
    errs = np.abs(errs) * 1e6
    check("'hmm' phase within 100 us on >= 5 of 6 records", np.sum(errs <= 100) >= 5,
          "errors (us): " + ", ".join(f"{e:.0f}" for e in errs))
    print(f"    SNR_max 3, vetoed  known phase: {summary(*f1s['known'])}\n"
          f"    SNR_max 3, vetoed  blind phase: {summary(*f1s['blind'])}")
    k, b = f1(*f1s["known"]), f1(*f1s["blind"])
    check("blind HMM-SF F1 within 0.03 of known-phase HMM-SF", b >= k - 0.03, f"known {k:.3f}, blind {b:.3f}")

    print()
    if _failures:
        print(f"FAILED ({len(_failures)}): " + "; ".join(_failures))
        return 1
    print("ALL CHECKS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
