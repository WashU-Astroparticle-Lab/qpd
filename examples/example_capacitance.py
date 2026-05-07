"""
Quantum-capacitance example.

Demonstrates the two new capabilities (issues #3 and #4):

  1. Compute and plot C_Q(n_g) for both parities.
  2. Generate a synthetic "measured" C_Q(n_g) trace with noise and a
     small charge offset, fit it, and overlay the recovered curve.
  3. Show the readout-basis chain
        measured Δφ  ->  Δf_r  ->  ΔC_Q
     for an off-resonant readout drive.
"""

import numpy as np
import matplotlib.pyplot as plt

from qpd import (
    QPD,
    phase_shift_to_frequency_shift,
    frequency_shift_to_quantum_capacitance,
)


def main():
    # ------------------------------------------------------------------
    # 1. Forward computation of C_Q(n_g)
    # ------------------------------------------------------------------
    e_j_hz = 8.335e9
    e_c_hz = 0.695e9        # E_J/E_C ≈ 12
    qpd = QPD(e_j_hz=e_j_hz, e_c_hz=e_c_hz)

    ng = np.linspace(0, 1, 401)
    fig1, _ = qpd.plot_quantum_capacitance(ng)
    fig1.savefig('cq_intrinsic.png', dpi=150, bbox_inches='tight')

    # ------------------------------------------------------------------
    # 2. Round-trip fit on a synthetic noisy trace
    # ------------------------------------------------------------------
    # Synthetic "measurement" in intrinsic units (Hz). Because we know
    # the data was generated in those units, we fix scale=1 — this
    # breaks the (E_C, scale) shape/amplitude degeneracy and lets the
    # fit pin down E_J and E_C tightly. With unknown absolute scale,
    # leave fit_scale=True; you'll then recover E_J/E_C and n_g0
    # robustly but only the product scale·E_C.
    ng0_true = 0.05
    _, cq_truth = qpd.compute_quantum_capacitance(ng - ng0_true)
    rng = np.random.default_rng(42)
    noise_amp = 0.02 * np.max(np.abs(cq_truth))
    cq_meas = cq_truth + rng.normal(0, noise_amp, size=ng.size)

    guess = QPD(e_j_hz=1.1 * e_j_hz, e_c_hz=0.9 * e_c_hz)
    fit = guess.fit_quantum_capacitance(
        ng, cq_meas, parity='odd', fit_scale=False, fixed_scale=1.0,
    )

    print('--- Fit result -----------------------------------------')
    print(f"  E_J  fitted = {fit['e_j_hz']/1e9:.4f} GHz   "
          f"(true {e_j_hz/1e9:.4f})")
    print(f"  E_C  fitted = {fit['e_c_hz']/1e9:.4f} GHz   "
          f"(true {e_c_hz/1e9:.4f})")
    print(f"  E_J/E_C     = {fit['ej_ec_ratio']:.3f}   "
          f"(true {e_j_hz/e_c_hz:.3f})")
    print(f"  n_g0 fitted = {fit['n_g0']:+.4f}        "
          f"(true {ng0_true:+.4f})")
    print(f"  scale       = {fit['scale']:.4g}")
    print('--------------------------------------------------------')

    fig2, _ = qpd.plot_capacitance_fit(ng, cq_meas, fit)
    fig2.savefig('cq_fit.png', dpi=150, bbox_inches='tight')

    # ------------------------------------------------------------------
    # 3. Readout chain: Δφ -> Δf_r -> ΔC_Q
    # ------------------------------------------------------------------
    f_r = 7.0e9          # bare resonator [Hz]
    kappa = 1.0e6        # 1 MHz linewidth
    c_r = 1.0e-13        # 0.1 pF effective resonator capacitance
    f_drive = f_r + 0.4 * kappa     # drive parked off-resonance

    # Suppose we measured a 50 mrad phase shift.
    delta_phi = 0.050
    detuning = f_drive - f_r
    df_r = phase_shift_to_frequency_shift(delta_phi, kappa, detuning)
    delta_cq = frequency_shift_to_quantum_capacitance(df_r, f_r, c_r)

    print('--- Readout chain --------------------------------------')
    print(f"  drive detuning Δ_d = {detuning/1e3:+.2f} kHz")
    print(f"  measured Δφ        = {delta_phi*1e3:.1f} mrad")
    print(f"  inferred Δf_r      = {df_r/1e3:+.3f} kHz")
    print(f"  inferred ΔC_Q      = {delta_cq:+.3e} F  "
          f"({delta_cq*1e18:+.3f} aF)")
    print('--------------------------------------------------------')

    plt.show()


if __name__ == '__main__':
    main()
