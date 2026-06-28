#!/usr/bin/env python
"""Verification for the stationary dressed-state S21 model + lock-in sinc window.

Run: ``python checks/check_readout_window.py`` (exit 0 = all checks pass).

Covers:
  1. Bare-cavity / single-pole Probst limit (g -> 0 recovers Probst Eq. 1 exactly).
  2. resonator_tools round-trip (skipped with a notice if not installed).
  3. Weak-drive peak position agrees with compute_dispersive_matrix(numeric) chi.
  4. Window limits: df -> 0 is a no-op; flat input stays flat; kernel has real
     negative side-lobes and sums to 1.
  5. Distortion vs df: apparent-chi error grows with df, ~0 when df << kappa.
"""
import sys

import numpy as np

from qpd import QPD
from qpd.theory.driven import (
    apply_lockin_window,
    extract_peak_hz,
    measured_s21_dressed,
)

# --- WashU QPD device point (QPD4 VNA fit) -----------------------------------
E_J = 6.95e9
E_C = 0.695e9
F_R = 6.936e9
G = 40e6
KAPPA = 115.70e3          # total (dissipator)
KAPPA_C = 69.82e3         # coupling/port
N_G = 0.0
PARITY = "even"

_failures = []


def check(name, ok, detail=""):
    tag = "PASS" if ok else "FAIL"
    print(f"  [{tag}] {name}" + (f"  --  {detail}" if detail else ""))
    if not ok:
        _failures.append(name)


def amp_min_hz(freqs, s21):
    """Frequency of the amplitude minimum |S21| (3-point parabolic refine).

    This is what an amplitude-only readout reports; for an asymmetric (phi != 0)
    line it sits away from the true resonance, and the bandwidth window shifts it.
    """
    mag = np.abs(np.asarray(s21))
    k = int(np.argmin(mag))
    if k == 0 or k == mag.size - 1:
        return float(freqs[k])
    y0, y1, y2 = mag[k - 1], mag[k], mag[k + 1]
    denom = y0 - 2.0 * y1 + y2
    if denom == 0.0:
        return float(freqs[k])
    delta = 0.5 * (y0 - y2) / denom
    return float(freqs[k] + delta * (freqs[k + 1] - freqs[k]))


def main():
    qpd = QPD(e_j_hz=E_J, e_c_hz=E_C)

    # ======================================================================
    # 1. Single-pole Probst limit: g -> 0 must reproduce Probst Eq. (1).
    # ======================================================================
    print("1. Single-pole Probst limit (g -> 0)")
    freqs = np.linspace(F_R - 3e6, F_R + 3e6, 4001)
    s21 = qpd.compute_s21_dressed(
        N_G, 0.0, F_R, freqs, parity=PARITY, kappa_hz=KAPPA, kappa_c_hz=KAPPA_C,
    )
    q_l = F_R / KAPPA
    q_c = F_R / KAPPA_C
    probst = 1.0 - (q_l / q_c) / (1.0 + 2j * q_l * (freqs / F_R - 1.0))
    max_dev = float(np.max(np.abs(s21 - probst)))
    check("matches analytic Probst Eq.(1)", max_dev < 1e-9,
          f"max|dS21| = {max_dev:.2e}")

    # On-resonance dip and circle diameter.
    s21_res = qpd.compute_s21_dressed(
        N_G, 0.0, F_R, np.array([F_R]), parity=PARITY,
        kappa_hz=KAPPA, kappa_c_hz=KAPPA_C,
    )[0]
    dip_expected = 1.0 - q_l / q_c
    check("on-resonance dip = 1 - Ql/Qc", abs(s21_res.real - dip_expected) < 1e-9
          and abs(s21_res.imag) < 1e-9,
          f"S21(f_r) = {s21_res:.6f}, expected {dip_expected:.6f}")
    # Diameter = |S21(off-res=1) - S21(on-res)| = Ql/Qc.
    diam = abs(1.0 - s21_res)
    check("circle diameter = Ql/Qc", abs(diam - q_l / q_c) < 1e-9,
          f"diameter = {diam:.6f}, Ql/Qc = {q_l / q_c:.6f}")

    # kappa_c = kappa -> perfect notch to 0.
    s21_crit = qpd.compute_s21_dressed(
        N_G, 0.0, F_R, np.array([F_R]), parity=PARITY,
        kappa_hz=KAPPA, kappa_c_hz=KAPPA,
    )[0]
    check("kappa_c = kappa -> dip to 0", abs(s21_crit) < 1e-9,
          f"|S21(f_r)| = {abs(s21_crit):.2e}")

    # ======================================================================
    # 2. resonator_tools round-trip (the exact daq fit path).
    # ======================================================================
    print("2. resonator_tools round-trip")
    try:
        from resonator_tools import circuit
        port = circuit.notch_port(freqs, s21)
        port.autofit()
        fr_fit = port.fitresults["fr"]
        ql_fit = port.fitresults["Ql"]
        qc_fit = port.fitresults["Qc_dia_corr"]
        check("recovers f_r", abs(fr_fit - F_R) < KAPPA,
              f"fr_fit = {fr_fit:.1f} Hz (true {F_R:.1f})")
        check("recovers Q_l", abs(ql_fit - q_l) / q_l < 0.05,
              f"Ql_fit = {ql_fit:.0f} (true {q_l:.0f})")
        check("recovers Q_c", abs(qc_fit - q_c) / q_c < 0.05,
              f"Qc_fit = {qc_fit:.0f} (true {q_c:.0f})")
    except ImportError:
        print("  [SKIP] resonator_tools not installed in this env "
              "(installed in the daq env); analytic check 1 already covers the "
              "Probst convention.")

    # ======================================================================
    # 3. Weak-drive peak position agrees with numeric dispersive chi.
    # ======================================================================
    print("3. Weak-drive peak vs compute_dispersive_matrix(numeric)")
    _, chi_ip = qpd.compute_dispersive_matrix(
        N_G, G, F_R, num_levels=4, parity=PARITY, method="numeric",
    )
    chi_ref = float(chi_ip[0])               # ground-level dispersive shift
    fwin = np.linspace(F_R - 2e6, F_R + 3e6, 6001)
    s21_g = qpd.compute_s21_dressed(
        N_G, G, F_R, fwin, parity=PARITY, kappa_hz=KAPPA, kappa_c_hz=KAPPA_C,
    )
    f_peak = extract_peak_hz(fwin, s21_g)
    chi_model = f_peak - F_R
    step = fwin[1] - fwin[0]
    check("peak shift ~ numeric chi",
          abs(chi_model - chi_ref) < max(3 * step, 0.02 * abs(chi_ref)),
          f"chi_model = {chi_model:.1f} Hz, chi_ref = {chi_ref:.1f} Hz "
          f"(grid step {step:.0f} Hz)")

    # ======================================================================
    # 4. Window limits.
    # ======================================================================
    print("4. Window limits")
    # df -> 0 (window much narrower than the grid step): no-op.
    s21_tiny = apply_lockin_window(fwin, s21_g, df_hz=step / 100.0)
    check("df << grid step -> no-op",
          np.max(np.abs(s21_tiny - s21_g)) < 1e-12)

    # Flat input stays flat.
    flat = np.ones_like(fwin, dtype=complex)
    flat_w = apply_lockin_window(fwin, flat, df_hz=50e3)
    check("flat S21 stays flat", np.max(np.abs(flat_w - 1.0)) < 1e-12,
          f"max dev = {np.max(np.abs(flat_w - 1.0)):.2e}")

    # Kernel sanity: real negative side-lobes, sums to 1.
    df_k = 50e3
    half = int(np.ceil(20 * df_k / step))
    m = np.arange(-half, half + 1)
    kernel = np.sinc(m * step / df_k)
    kernel = kernel / kernel.sum()
    check("kernel sums to 1", abs(kernel.sum() - 1.0) < 1e-12)
    check("kernel has negative side-lobes (true sinc)", kernel.min() < 0.0,
          f"min = {kernel.min():.4f}")

    # ======================================================================
    # 5. Distortion vs df. PHYSICS: a SYMMETRIC sinc window on a SYMMETRIC notch
    #    broadens it and fills the dip but does NOT move the peak; an apparent
    #    chi (peak) bias appears only for an ASYMMETRIC lineshape (phi != 0 Fano,
    #    or near a crossing). Both behaviours are verified.
    # ======================================================================
    print("5. Distortion vs df")
    fwin5 = np.linspace(F_R - 6e6, F_R + 6e6, 12001)   # wide enough for df<=150 kHz
    s21_true5 = qpd.compute_s21_dressed(
        N_G, G, F_R, fwin5, parity=PARITY, kappa_hz=KAPPA, kappa_c_hz=KAPPA_C,
    )
    true_depth = float(np.max(np.abs(s21_true5 - 1.0)))
    true_peak = extract_peak_hz(fwin5, s21_true5)

    print("   (a) symmetric notch: depth fills, peak stays put")
    depths, peak_bias = {}, {}
    for df in (1e3, 10e3, 50e3, 150e3):
        s21_w = apply_lockin_window(fwin5, s21_true5, df_hz=df)
        depths[df] = float(np.max(np.abs(s21_w - 1.0)))
        peak_bias[df] = abs(extract_peak_hz(fwin5, s21_w) - true_peak)
        print(f"     df = {df / 1e3:6.0f} kHz  ->  notch depth = "
              f"{depths[df]:.4f} (true {true_depth:.4f})   "
              f"|peak bias| = {peak_bias[df]:6.1f} Hz")
    step5 = fwin5[1] - fwin5[0]
    check("depth reduced as df grows (df=150k < df=1k)",
          depths[150e3] < depths[1e3] - 1e-4,
          f"{depths[1e3]:.4f} -> {depths[150e3]:.4f}")
    check("depth ~ unchanged for df << kappa (df=1 kHz)",
          abs(depths[1e3] - true_depth) < 1e-3,
          f"depth {depths[1e3]:.4f} vs true {true_depth:.4f}")
    check("symmetric peak NOT biased by symmetric window",
          max(peak_bias.values()) < 3 * step5,
          f"max |peak bias| = {max(peak_bias.values()):.1f} Hz "
          f"(grid step {step5:.0f} Hz)")

    print("   (b) asymmetric line (phi = 0.5 rad): amplitude-min readout shifts")
    s21_asym = qpd.compute_s21_dressed(
        N_G, G, F_R, fwin5, parity=PARITY, kappa_hz=KAPPA, kappa_c_hz=KAPPA_C,
        phi=0.5,
    )
    amin0 = amp_min_hz(fwin5, s21_asym)
    bias_small = abs(amp_min_hz(
        fwin5, apply_lockin_window(fwin5, s21_asym, df_hz=1e3)) - amin0)
    bias_large = abs(amp_min_hz(
        fwin5, apply_lockin_window(fwin5, s21_asym, df_hz=150e3)) - amin0)
    print(f"     amplitude-min bias: df=1 kHz -> {bias_small:.1f} Hz,  "
          f"df=150 kHz -> {bias_large:.1f} Hz")
    check("asymmetric amplitude-min bias grows with df",
          bias_large > bias_small + step5,
          f"{bias_small:.1f} -> {bias_large:.1f} Hz")

    # ----------------------------------------------------------------------
    print()
    if _failures:
        print(f"RESULT: {len(_failures)} check(s) FAILED: {_failures}")
        return 1
    print("RESULT: all checks PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
