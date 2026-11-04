#!/usr/bin/env python
"""Verify the LED flip network, ``qpd.reconstruction.led_network.reconstruct_parity_flips_led``.

Covers the API contract (input forms, the t0 axis, the post-onset veto, the derived branch, the estimated LED phase
when no trigger time is given, errors and the HMM-SF fallback), the bundled weights, accuracy on records from the
LED-pulsing simulator against the HMM-SF baseline, the high-SNR guard, the coloured-noise guard (auto-whitening,
and the hard flag when whitening is switched off), all of these again for 10 kHz records, LED schedules other than the
trained one (0.1 ms on / 60 ms off), and other sampling rates.

Needs the optional extra (``pip install "qpd[led]"``); without torch it prints SKIP and exits 0.
Run: ``python checks/check_led_network.py`` (exit 0 = all checks pass). Seeds are fixed.
"""

import hashlib
import os
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

from check_led_hmm import FS, PERIOD, RES, TOL, FastSimulator, max_separation, simulate  # noqa: E402
from qpd import QPD  # noqa: E402
from qpd.simulator import LEDPulseModel, WhiteGaussianNoise  # noqa: E402
from qpd.mlebench.grade import match_times  # noqa: E402
from qpd.reconstruction.led_network import LEDReconstructionResult, reconstruct_parity_flips_hmmsf, reconstruct_parity_flips_led  # noqa: E402
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


def simulate_schedule(seed, snr_max, on_time, off_time, duration=2.0, n_qp=3.0, rate=FS):
    """Like check_led_hmm.simulate, on another LED schedule or (at the same noise density) a faster sampling rate."""
    first = float(np.random.default_rng(seed).uniform(0, on_time + off_time))
    sim = FastSimulator(
        qpd=QPD(e_j_hz=6.95e9, e_c_hz=0.695e9, coupling_g_hz=40e6), resonator=RES,
        f_drive=RES.f_r, sample_rate=rate, noise=WhiteGaussianNoise(sigma=max_separation() / snr_max * np.sqrt(rate / FS)),
        led_pulsing=LEDPulseModel(on_time=on_time, off_time=off_time, first_pulse_time=first, expected_n_qp=n_qp, charge_jump_rate_hz=0.0))
    res = sim.simulate(duration=duration, seed=seed)
    return res.iq, res.flip_times, first


def slower(iq, rate):
    """The same record sampled at ``rate`` < 100 kHz: every sample averages its full sample period."""
    n = int(iq.size * rate / FS)
    k = np.floor(np.arange(iq.size) * rate / FS).astype(np.int64)
    keep = k < n
    return (np.bincount(k[keep], iq.real[keep], n) + 1j * np.bincount(k[keep], iq.imag[keep], n)) / np.bincount(k[keep], minlength=n)


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
    check("no schedule flag on the trained schedule", "untrained_schedule" not in r.flags and "untested_schedule" not in r.flags)
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
    raises("a non-positive sample rate raises", lambda: reconstruct_parity_flips_led(iq, 0.0, first_pulse_time_s=fp))

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
          and r3.diagnostics["resampled"]["n_samples"] == iq10.size * 10)
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

    print("LED 0.1 ms on / 60 ms off")
    ON, OFF = 1e-4, 60e-3
    SCHEDULE = dict(on_time_s=ON, off_time_s=OFF)
    from qpd.mlebench.led_hmm import predict_led_hmm as hmm_sf
    from qpd.reconstruction.led_network.health import onset_times, veto_mask
    SF_KW = dict(on_time_s=ON, off_time_s=OFF, pin_quiet_rate=True, smooth_fit_s=2e-3, min_smoothed_contrast=2.5)
    for rate, label in ((1e4, "10 kHz"), (FS, "100 kHz")):
        f1_net, f1_sf = [], []
        for seed, snr, n_qp in [(31, 3.0, 3.0), (32, 1.0, 3.0), (33, 2.0, 10.0), (34, 1.0, 1.0)]:
            iq, truth, fp = simulate_schedule(seed, snr, ON, OFF, n_qp=n_qp)
            r = reconstruct_parity_flips_led(to_10khz(iq) if rate == 1e4 else iq, rate, first_pulse_time_s=fp, **SCHEDULE)
            e, p, f = eff_pur(truth, r.flip_times)
            _, _, fs = eff_pur(truth, r.hmmsf_flip_times)
            f1_net.append(f); f1_sf.append(fs)
            check(f"{label}, seed {seed} (SNR {snr}, {n_qp:g} per pulse): network, no hard flag, efficiency >= 0.5, purity >= 0.75",
                  r.decoder == "network" and not r.degenerate and e >= 0.5 and p >= 0.75, f"eff {e:.3f} pur {p:.3f} F1 {f:.3f} (HMM-SF F1 {fs:.3f})")
        check(f"{label}: mean F1, network >= HMM-SF", np.mean(f1_net) >= np.mean(f1_sf), f"{np.mean(f1_net):.3f} vs {np.mean(f1_sf):.3f}")
    iq, truth, fp = simulate_schedule(35, 3.0, ON, OFF, duration=1.0)
    iq10 = to_10khz(iq)
    r = reconstruct_parity_flips_led(iq10, 1e4, first_pulse_time_s=fp, **SCHEDULE)
    check("soft flag untrained_schedule fires, reported on the 10 kHz grid", "untrained_schedule" in r.flags and "untrained_schedule" not in r.diagnostics["hard_flags_fired"]
          and r.sample_rate == 1e4 and r.n_samples == iq10.size, str(r.flags))
    ons = onset_times(fp, iq10.size / 1e4, ON + OFF)
    check("the veto follows the given onsets", ons.size >= 16 and np.all(veto_mask(r.flip_times, ons)))
    native, _ = hmm_sf(iq10, 1e4, first_pulse_time_s=fp, **SF_KW)
    check("HMM-SF baseline uses the given schedule", np.array_equal(r.hmmsf_flip_times, np.sort(native[veto_mask(native, ons)])))
    rb = reconstruct_parity_flips_led(iq10, 1e4, **SCHEDULE)
    err = ((rb.first_pulse_time_s - fp + (ON + OFF) / 2) % (ON + OFF)) - (ON + OFF) / 2
    check("first_pulse_time_s omitted: phase estimated within 0.3 ms", rb.phase_estimated and abs(err) < 3e-4, f"error {1e6*err:+.0f} us")
    iq, truth, fp = simulate_schedule(36, 20.0, ON, OFF, duration=1.0)
    rf = reconstruct_parity_flips_led(iq, FS, first_pulse_time_s=fp, fallback_hmmsf=True, **SCHEDULE)
    single = iq.real.astype(np.float32).astype(float) + 1j * iq.imag.astype(np.float32).astype(float)     # as preprocess casts it
    native, _ = hmm_sf(single, FS, first_pulse_time_s=fp, **SF_KW)
    native = np.sort(native[veto_mask(native, onset_times(fp, iq.size / FS, ON + OFF))])
    _, p, _ = eff_pur(truth, rf.flip_times)
    check("100 kHz fallback runs HMM-SF on the given schedule", rf.decoder == "hmm_sf_fallback" and np.array_equal(rf.flip_times, native) and p >= 0.9, f"purity {p:.3f}")
    iq, truth, fp = simulate_schedule(37, 3.0, ON, 20e-3, duration=1.0)
    ru = reconstruct_parity_flips_led(iq, FS, first_pulse_time_s=fp, on_time_s=ON, off_time_s=20e-3)
    check("a period below the tested range raises the hard flag", ru.degenerate and "untested_schedule" in ru.flags)
    raises("a non-positive on time raises", lambda: reconstruct_parity_flips_led(iq, FS, first_pulse_time_s=fp, on_time_s=0.0))

    print("other sampling rates")
    from qpd.reconstruction.led_network.inputs import average_100khz, rebuild_100khz
    iq, truth, fp = simulate_schedule(41, 3.0, ON, OFF, duration=1.0, rate=2.5e5)
    y_i, _ = average_100khz(iq.real, iq.imag, 2.5e5)
    check("averaging 250 kHz gives 100 kHz samples with the right mean", y_i.size == 100000 and np.isclose(y_i.mean(), iq.real.mean(), atol=1e-3))
    r = reconstruct_parity_flips_led(iq, 2.5e5, first_pulse_time_s=fp, **SCHEDULE)
    e, p, f = eff_pur(truth, r.flip_times)
    check("250 kHz: network, reported on the record's grid, efficiency >= 0.6, purity >= 0.8", r.decoder == "network" and not r.degenerate and r.sample_rate == 2.5e5
          and r.n_samples == iq.size and r.branch.shape == iq.shape and e >= 0.6 and p >= 0.8, f"eff {e:.3f} pur {p:.3f}")
    iq, truth, fp = simulate_schedule(42, 3.0, ON, OFF, duration=1.0)
    for rate in (5e4, 3e4, 5e3):
        x = slower(iq, rate)
        y_i, _ = rebuild_100khz(x.real, x.imag, rate)
        edges = np.round(np.arange(x.size + 1) * FS / rate).astype(np.int64)
        check(f"{rate/1e3:g} kHz: rebuilt blocks keep each sample's mean", y_i.size == edges[-1] and np.allclose(np.add.reduceat(y_i, edges[:-1]) / np.diff(edges), x.real))
        r = reconstruct_parity_flips_led(x, rate, first_pulse_time_s=fp, **SCHEDULE)
        e, p, f = eff_pur(truth, r.flip_times)
        check(f"{rate/1e3:g} kHz: network, reported on the record's grid, efficiency >= 0.6, purity >= 0.8", r.decoder == "network" and not r.degenerate and r.sample_rate == rate
              and r.n_samples == x.size and r.branch.shape == x.shape and e >= 0.6 and p >= 0.8, f"eff {e:.3f} pur {p:.3f}")
    rb = reconstruct_parity_flips_led(slower(iq, 3e4), 3e4, **SCHEDULE)
    err = ((rb.first_pulse_time_s - fp + (ON + OFF) / 2) % (ON + OFF)) - (ON + OFF) / 2
    check("30 kHz, first_pulse_time_s omitted: phase estimated within 0.3 ms", rb.phase_estimated and abs(err) < 3e-4, f"error {1e6*err:+.0f} us")
    r = reconstruct_parity_flips_led(slower(iq, 1e3), 1e3, first_pulse_time_s=fp, **SCHEDULE)
    check("1 kHz runs and reports on the record's grid", r.sample_rate == 1e3 and r.branch.shape == (slower(iq, 1e3).size,), f"{r.flip_times.size} flips of {truth.size}")

    print("HMM-SF on its own (reconstruct_parity_flips_hmmsf)")
    iq, truth, fp = simulate_schedule(43, 3.0, ON, OFF, duration=1.0)
    fast, fast_truth, fast_fp = simulate_schedule(44, 3.0, ON, OFF, duration=1.0, rate=2.5e5)
    for label, x, rate, first in (("100 kHz", iq, FS, fp), ("10 kHz", slower(iq, 1e4), 1e4, fp), ("30 kHz", slower(iq, 3e4), 3e4, fp), ("250 kHz", fast, 2.5e5, fast_fp)):
        h = reconstruct_parity_flips_hmmsf(x, rate, first, **SCHEDULE)
        n = reconstruct_parity_flips_led(x, rate, first, **SCHEDULE)
        hb = reconstruct_parity_flips_hmmsf(x, rate, **SCHEDULE)
        nb = reconstruct_parity_flips_led(x, rate, **SCHEDULE)
        check(f"{label}: same flips as the network call's hmmsf_flip_times, trigger time given and estimated", np.array_equal(h.flip_times, n.hmmsf_flip_times)
              and np.array_equal(hb.flip_times, nb.hmmsf_flip_times) and hb.first_pulse_time_s == nb.first_pulse_time_s and hb.phase_estimated and not h.phase_estimated)
        check(f"{label}: reported on the record's grid, decoder hmm_sf, no confidence", h.decoder == "hmm_sf" and h.sample_rate == rate and h.n_samples == x.size
              and h.branch.shape == x.shape and np.isnan(h.confidence).all() and h.confidence.size == h.flip_times.size and not h.degenerate)
    h = reconstruct_parity_flips_hmmsf(iq, FS, fp, **SCHEDULE)
    e, p, f = eff_pur(truth, h.flip_times)
    check("efficiency >= 0.6, purity >= 0.8", e >= 0.6 and p >= 0.8, f"eff {e:.3f} pur {p:.3f} F1 {f:.3f}")
    d = reconstruct_parity_flips_hmmsf(simulate(21, 3.0, duration=1.0, n_qp=3.0)[0], FS, simulate(21, 3.0, duration=1.0, n_qp=3.0)[2])
    r = reconstruct_parity_flips_led(simulate(21, 3.0, duration=1.0, n_qp=3.0)[0], FS, simulate(21, 3.0, duration=1.0, n_qp=3.0)[2])
    check("schedule defaults to the trained one", np.array_equal(d.flip_times, r.hmmsf_flip_times))
    s = reconstruct_parity_flips_hmmsf(iq, FS, fp + 10.0, t0=10.0, **SCHEDULE)
    check("t0 shifts the flip times and nothing else", np.allclose(s.flip_times - 10.0, h.flip_times, rtol=0, atol=1e-9))
    rng = np.random.default_rng(43)
    x = slower(iq, 1e4)
    e = (rng.normal(size=x.size) + 1j * rng.normal(size=x.size)) * np.std(np.diff(x.real)) / np.sqrt(2)
    coloured = x + lfilter([np.sqrt(1 - 0.3 ** 2)], [1, -0.3], e)
    ca, cn = reconstruct_parity_flips_hmmsf(coloured, 1e4, fp, **SCHEDULE), reconstruct_parity_flips_led(coloured, 1e4, fp, **SCHEDULE)
    check("coloured noise: whitened, same flips as the network call", ca.diagnostics["whitening"]["applied"] and np.array_equal(ca.flip_times, cn.hmmsf_flip_times))
    co = reconstruct_parity_flips_hmmsf(coloured, 1e4, fp, whiten="off", **SCHEDULE)
    check("whiten='off' on coloured noise raises the hard flag", co.degenerate and "coloured_noise" in co.flags)
    rng = np.random.default_rng(45)
    noise = rng.normal(size=20000) + 1j * rng.normal(size=20000)
    nn = reconstruct_parity_flips_hmmsf(noise, 1e4, 0.01, **SCHEDULE)
    check("pure noise: flips only if a window passes the gate, hard flag if none does", nn.degenerate == (nn.diagnostics["hmmsf"]["info"]["dead_segments"] >= nn.diagnostics["hmmsf"]["info"]["segments"]),
          f"{nn.flip_times.size} flips, degenerate {nn.degenerate}")
    raises("a non-positive sample rate raises", lambda: reconstruct_parity_flips_hmmsf(iq, 0.0, fp))
    raises("a non-positive on time raises", lambda: reconstruct_parity_flips_hmmsf(iq, FS, fp, on_time_s=0.0))
    import subprocess
    code = ("import sys\n"
            "class NoTorch:\n"
            "    def find_spec(self, name, path=None, target=None):\n"
            "        if name.split('.')[0] == 'torch':\n"
            "            raise ImportError(name)\n"
            "sys.meta_path.insert(0, NoTorch())\n"
            "import numpy as np\n"
            "from qpd.reconstruction import reconstruct_parity_flips_hmmsf, reconstruct_parity_flips_led\n"
            "rng = np.random.default_rng(0); iq = rng.normal(size=20000) + 1j * rng.normal(size=20000)\n"
            "print(reconstruct_parity_flips_hmmsf(iq, 1e4, 0.01, on_time_s=1e-4, off_time_s=60e-3).decoder)\n"
            "try:\n    reconstruct_parity_flips_led(iq, 1e4, 0.01, on_time_s=1e-4, off_time_s=60e-3)\n"
            "except ImportError as error:\n    print('ImportError', 'pip install' in str(error))\n")
    src = str(Path(__file__).resolve().parents[1] / "src")
    p = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, env={**os.environ, "PYTHONPATH": src})
    check("without torch: HMM-SF runs, the network raises a clear ImportError", p.stdout.split() == ["hmm_sf", "ImportError", "True"], p.stdout.strip() + p.stderr[-200:])

    print()
    if _failures:
        print(f"FAILED {len(_failures)} check(s): {_failures}")
        sys.exit(1)
    print("All LED-network checks passed.")


if __name__ == "__main__":
    main()
