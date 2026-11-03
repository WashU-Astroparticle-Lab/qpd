#!/usr/bin/env python
"""Verify the LED flip network, ``qpd.reconstruction.led_network.reconstruct_parity_flips_led``.

Covers the API contract (input forms, the t0 axis, the post-onset veto, the derived branch, the estimated LED phase
when no trigger time is given, errors and the HMM-SF fallback), the bundled weights, accuracy on records from the
LED-pulsing simulator against the HMM-SF baseline, the high-SNR guard, the coloured-noise guard (auto-whitening,
and the hard flag when whitening is switched off), and all of these again for 10 kHz records.

Needs the optional extra (``pip install "qpd[led]"``); without torch it prints SKIP and exits 0.
Run: ``python checks/check_led_network.py`` (exit 0 = all checks pass). Seeds are fixed.
"""

import hashlib
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

try:
    import torch  # noqa: F401
except ImportError:
    print("SKIP: torch is not installed (pip install \"qpd[led]\")")
    sys.exit(0)

from check_led_hmm import FS, PERIOD, TOL, simulate  # noqa: E402
from qpd.mlebench.grade import match_times  # noqa: E402
from qpd.reconstruction.led_network import LEDReconstructionResult, reconstruct_parity_flips_led  # noqa: E402
from qpd.reconstruction.led_network.reconstruct import DEFAULT_WEIGHTS  # noqa: E402

WEIGHTS_SHA256 = "928b568b6e2817f3ea1e3cecd5f099cfa6c7b815b41453e3f30ab7f639abe63e"
VETO = 300e-6
_failures = []


def check(name, ok, detail=""):
    tag = "PASS" if ok else "FAIL"
    print(f"  [{tag}] {name}" + (f"  --  {detail}" if detail else ""))
    if not ok:
        _failures.append(name)


def raises(name, action, exc=ValueError):
    try:
        action()
    except exc:
        check(name, True)
    else:
        check(name, False)


def eff_pur(truth, pred):
    m = match_times(np.sort(truth), np.sort(pred), TOL)
    eff, pur = m.tp / max(len(truth), 1), m.tp / max(len(pred), 1)
    return eff, pur, 2 * eff * pur / (eff + pur) if eff + pur else 0.0


def to_10khz(iq):
    """The same record sampled at 10 kHz: every sample averages its full 100 us."""
    return iq[:iq.size // 10 * 10].reshape(-1, 10).mean(axis=1)


def main():
    print("weights")
    sha = hashlib.sha256(DEFAULT_WEIGHTS.read_bytes()).hexdigest()
    check("bundled weights file is the released one", sha == WEIGHTS_SHA256, sha[:12])
    from qpd.reconstruction.led_network.network import N_PARAMETERS, load_network
    model = load_network(DEFAULT_WEIGHTS)
    check("network holds 9848 learned values", sum(p.numel() for p in model.parameters()) == N_PARAMETERS == 9848)

    print("API contract")
    iq, truth, fp = simulate(21, 3.0, duration=1.0, n_qp=3.0)
    r = reconstruct_parity_flips_led(iq, FS, first_pulse_time_s=fp)
    check("returns LEDReconstructionResult", isinstance(r, LEDReconstructionResult))
    check("flip times sorted and inside the record", np.all(np.diff(r.flip_times) >= 0)
          and np.all((r.flip_times >= 0) & (r.flip_times < iq.size / FS)))
    check("confidence per flip, in [0.5, 1]", r.confidence.shape == r.flip_times.shape
          and np.all((r.confidence >= 0.5) & (r.confidence <= 1)))
    onsets = np.arange(fp, iq.size / FS, PERIOD)
    k = np.searchsorted(onsets, r.flip_times, side="right") - 1
    since = np.where(k >= 0, r.flip_times - onsets[np.maximum(k, 0)], np.inf)
    check("no flip within 300 us after an LED onset", np.all(since >= VETO))
    check("decoder is the network on a healthy record", r.decoder == "network" and not r.degenerate, str(r.flags))
    rq = reconstruct_parity_flips_led((iq.real, iq.imag), FS, first_pulse_time_s=fp)
    check("(I, Q) input equals complex input", np.array_equal(rq.flip_times, r.flip_times))
    rs = reconstruct_parity_flips_led(iq, FS, first_pulse_time_s=fp + 10.0, t0=10.0)
    check("t0 shifts the flip times and nothing else", np.allclose(rs.flip_times - 10.0, r.flip_times, rtol=0, atol=1e-9)
          and rs.flip_times.size == r.flip_times.size)
    b = r.branch
    check("branch toggles once per flip", b.shape == (iq.size,) and b.dtype == np.int8
          and int(np.count_nonzero(np.diff(b))) == r.flip_times.size - int(np.sum(np.round(r.flip_times * FS) == 0)))
    check("HMM-SF baseline flips reported", r.hmmsf_flip_times is not None and r.hmmsf_flip_times.size > 0)
    check("rate_hz and mu_hat are finite", np.isfinite(r.rate_hz) and np.isfinite(r.mu_hat))
    check("known trigger time is reported back", not r.phase_estimated and r.first_pulse_time_s == fp)
    rb = reconstruct_parity_flips_led(iq, FS)
    err = ((rb.first_pulse_time_s - fp + PERIOD / 2) % PERIOD) - PERIOD / 2
    _, _, f_known = eff_pur(truth, r.flip_times)
    _, _, f_blind = eff_pur(truth, rb.flip_times)
    check("first_pulse_time_s omitted: phase estimated within 0.3 ms", rb.phase_estimated and abs(err) < 3e-4, f"error {1e6*err:+.0f} us")
    check("estimated-phase decode within 0.05 F1 of the known-phase decode", f_blind >= f_known - 0.05,
          f"F1 {f_blind:.3f} vs {f_known:.3f}")
    rbs = reconstruct_parity_flips_led(iq, FS, t0=10.0)
    check("estimated phase reported on the t0 axis", np.isclose(rbs.first_pulse_time_s - 10.0, rb.first_pulse_time_s, rtol=0, atol=1e-9)
          and np.allclose(rbs.flip_times - 10.0, rb.flip_times, rtol=0, atol=1e-9))
    raises("unsupported sample rate raises", lambda: reconstruct_parity_flips_led(iq, 2e5, first_pulse_time_s=fp))
    rf = reconstruct_parity_flips_led(iq, 2e5, first_pulse_time_s=fp, fallback_hmmsf=True)
    check("unsupported rate with fallback_hmmsf uses HMM-SF", rf.decoder == "hmm_sf_fallback" and rf.degenerate)

    print("accuracy on simulated LED records (network vs HMM-SF)")
    f1_net, f1_sf = [], []
    for seed, snr, n_qp in [(11, 3.0, 3.0), (12, 1.5, 3.0), (13, 3.0, 10.0), (14, 1.0, 1.0)]:
        iq, truth, fp = simulate(seed, snr, duration=2.0, n_qp=n_qp)
        r = reconstruct_parity_flips_led(iq, FS, first_pulse_time_s=fp)
        e, p, f = eff_pur(truth, r.flip_times)
        _, _, fs = eff_pur(truth, r.hmmsf_flip_times)
        f1_net.append(f); f1_sf.append(fs)
        check(f"seed {seed} (SNR {snr}, {n_qp:g} per pulse): no hard flag, efficiency >= 0.6, purity >= 0.85",
              not r.degenerate and e >= 0.6 and p >= 0.85, f"eff {e:.3f} pur {p:.3f} F1 {f:.3f} (HMM-SF F1 {fs:.3f})")
    check("mean F1: network >= HMM-SF", np.mean(f1_net) >= np.mean(f1_sf), f"{np.mean(f1_net):.3f} vs {np.mean(f1_sf):.3f}")

    print("above the trained SNR range")
    iq, truth, fp = simulate(16, 20.0, duration=1.0, n_qp=10.0)
    rh = reconstruct_parity_flips_led(iq, FS, first_pulse_time_s=fp)
    check("SNR 20 record is flagged degenerate", rh.degenerate and "snr_above_trained_range" in rh.flags,
          f"snr_hat {rh.diagnostics['stats']['snr_hat']:.1f}")
    rf = reconstruct_parity_flips_led(iq, FS, first_pulse_time_s=fp, fallback_hmmsf=True)
    _, _, f_sf = eff_pur(truth, rf.flip_times)
    check("fallback_hmmsf returns the HMM-SF flips there", rf.decoder == "hmm_sf_fallback"
          and np.array_equal(rf.flip_times, rh.hmmsf_flip_times), f"HMM-SF F1 {f_sf:.3f}")

    print("coloured noise")
    from scipy.signal import lfilter
    iq, truth, fp = simulate(15, 3.0, duration=1.0, n_qp=3.0)
    rng = np.random.default_rng(15)
    e = (rng.normal(size=iq.size) + 1j * rng.normal(size=iq.size)) * np.std(np.diff(iq.real)) / np.sqrt(2)
    coloured = iq + lfilter([np.sqrt(1 - 0.3 ** 2)], [1, -0.3], e)
    ra = reconstruct_parity_flips_led(coloured, FS, first_pulse_time_s=fp)
    check("auto-whitening applied for AR(1) noise", ra.diagnostics["whitening"]["applied"],
          f"rho_hat {ra.diagnostics['whitening']['rho_hat_raw']:.3f}")
    check("no coloured-noise flag after whitening", "coloured_noise" not in ra.flags)
    ro = reconstruct_parity_flips_led(coloured, FS, first_pulse_time_s=fp, whiten="off")
    check("whiten='off' on coloured noise raises the hard flag", ro.degenerate and "coloured_noise" in ro.flags)

    print("10 kHz records")
    iq, truth, fp = simulate(21, 3.0, duration=1.0, n_qp=3.0)
    iq10 = to_10khz(iq)
    r = reconstruct_parity_flips_led(iq10, 1e4, first_pulse_time_s=fp)
    check("decoded by the network, reported on the 10 kHz grid", r.decoder == "network" and not r.degenerate
          and r.sample_rate == 1e4 and r.n_samples == iq10.size and r.branch.shape == iq10.shape, str(r.flags))
    check("flip times sorted and inside the record", np.all(np.diff(r.flip_times) >= 0)
          and np.all((r.flip_times >= 0) & (r.flip_times < iq10.size / 1e4)))
    onsets = np.arange(fp, iq10.size / 1e4, PERIOD)
    k = np.searchsorted(onsets, r.flip_times, side="right") - 1
    since = np.where(k >= 0, r.flip_times - onsets[np.maximum(k, 0)], np.inf)
    check("no flip within 300 us after an LED onset", np.all(since >= VETO))
    r2 = reconstruct_parity_flips_led((iq10.real, iq10.imag), 1e4, first_pulse_time_s=fp)
    check("same rebuild_seed reproduces the output", np.array_equal(r2.flip_times, r.flip_times))
    r3 = reconstruct_parity_flips_led(iq10, 1e4, first_pulse_time_s=fp, rebuild_seed=1)
    check("the block noise follows rebuild_seed", r3.diagnostics["options"]["rebuild_seed"] == 1
          and r3.diagnostics["rebuild"]["n_samples"] == iq10.size * 10)
    rs = reconstruct_parity_flips_led(iq10, 1e4, first_pulse_time_s=fp + 10.0, t0=10.0)
    check("t0 shifts the flip times and nothing else", np.allclose(rs.flip_times - 10.0, r.flip_times, rtol=0, atol=1e-9)
          and np.allclose(rs.hmmsf_flip_times - 10.0, r.hmmsf_flip_times, rtol=0, atol=1e-9))
    from qpd.mlebench.led_hmm import predict_led_hmm
    native, _ = predict_led_hmm(iq10, 1e4, first_pulse_time_s=fp, off_time_s=PERIOD - 50e-6, pin_quiet_rate=True,
                                smooth_fit_s=2e-3, min_smoothed_contrast=2.5)
    since = (native - fp) % PERIOD
    native = np.sort(native[(native < fp) | (since >= VETO)])
    check("HMM-SF baseline runs on the 10 kHz samples themselves", np.array_equal(r.hmmsf_flip_times, native))
    rb = reconstruct_parity_flips_led(iq10, 1e4)
    err = ((rb.first_pulse_time_s - fp + PERIOD / 2) % PERIOD) - PERIOD / 2
    check("first_pulse_time_s omitted: phase estimated within 0.3 ms", rb.phase_estimated and abs(err) < 3e-4, f"error {1e6*err:+.0f} us")

    f1_net, f1_sf = [], []
    for seed, snr, n_qp in [(11, 3.0, 3.0), (12, 1.5, 3.0), (13, 3.0, 10.0), (14, 1.0, 1.0)]:
        iq, truth, fp = simulate(seed, snr, duration=2.0, n_qp=n_qp)
        r = reconstruct_parity_flips_led(to_10khz(iq), 1e4, first_pulse_time_s=fp)
        e, p, f = eff_pur(truth, r.flip_times)
        _, _, fs = eff_pur(truth, r.hmmsf_flip_times)
        f1_net.append(f); f1_sf.append(fs)
        check(f"seed {seed} (100 kHz-equivalent SNR {snr}, {n_qp:g} per pulse): no hard flag, efficiency >= 0.55, purity >= 0.8",
              not r.degenerate and e >= 0.55 and p >= 0.8, f"eff {e:.3f} pur {p:.3f} F1 {f:.3f} (HMM-SF F1 {fs:.3f})")
    check("mean F1: network >= HMM-SF at 10 kHz", np.mean(f1_net) >= np.mean(f1_sf), f"{np.mean(f1_net):.3f} vs {np.mean(f1_sf):.3f}")

    iq, truth, fp = simulate(16, 20.0, duration=1.0, n_qp=10.0)
    iq10 = to_10khz(iq)
    rh = reconstruct_parity_flips_led(iq10, 1e4, first_pulse_time_s=fp, fallback_hmmsf=True)
    check("above the trained SNR: flagged, HMM-SF answers on the 10 kHz samples", rh.degenerate and "snr_above_trained_range" in rh.flags
          and rh.decoder == "hmm_sf_fallback" and np.array_equal(rh.flip_times, rh.hmmsf_flip_times), str(rh.flags))

    iq, truth, fp = simulate(15, 3.0, duration=1.0, n_qp=3.0)
    iq10 = to_10khz(iq)
    rng = np.random.default_rng(15)
    e = (rng.normal(size=iq10.size) + 1j * rng.normal(size=iq10.size)) * np.std(np.diff(iq10.real)) / np.sqrt(2)
    coloured = iq10 + lfilter([np.sqrt(1 - 0.3 ** 2)], [1, -0.3], e)
    ra = reconstruct_parity_flips_led(coloured, 1e4, first_pulse_time_s=fp)
    check("coloured 10 kHz noise: measured at 10 kHz and whitened", ra.diagnostics["whitening"]["applied"] and "coloured_noise" not in ra.flags,
          f"rho_hat {ra.diagnostics['whitening']['rho_hat_raw']:.3f}")
    ro = reconstruct_parity_flips_led(coloured, 1e4, first_pulse_time_s=fp, whiten="off")
    check("whiten='off' on coloured 10 kHz noise raises the hard flag", ro.degenerate and "coloured_noise" in ro.flags)
    rw = reconstruct_parity_flips_led(iq10, 1e4, first_pulse_time_s=fp)
    check("white 10 kHz noise is left unwhitened", not rw.diagnostics["whitening"]["applied"] and "coloured_noise" not in rw.flags,
          f"rho_hat {rw.diagnostics['whitening']['rho_hat_raw']:.3f}")

    print()
    if _failures:
        print(f"FAILED {len(_failures)} check(s): {_failures}")
        sys.exit(1)
    print("All LED-network checks passed.")


if __name__ == "__main__":
    main()
