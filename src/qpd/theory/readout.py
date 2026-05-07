"""
Resonator-basis readout helpers.

Pure utility functions that map between three measured / inferred quantities:

    measured phase shift  <-->  resonator frequency shift  <-->  ΔC_Q

These are independent of any QPD instance state and are intended to be used
together with `QPD.compute_quantum_capacitance` to interpret experimental
phase / frequency-shift data.

Conventions
-----------
- Frequencies are in Hz, phases in radians, capacitances in Farads.
- A driven Lorentzian cavity has phase response

      φ(Δ_d) = arctan(2 Δ_d / κ),     Δ_d = f_drive - f_r,

  where κ is the loaded resonator linewidth (FWHM, in Hz). When the qubit
  shifts the bare resonator frequency by Δf (i.e. f_r → f_r + Δf), the
  drive detuning becomes Δ_d - Δf and the measured phase becomes

      φ_new = arctan(2 (Δ_d - Δf) / κ).

  `phase_shift_to_frequency_shift` inverts this in closed form for a
  measured *change* in phase Δφ = φ_new - φ_old:

      Δf = Δ_d - (κ/2) · tan(arctan(2 Δ_d / κ) + Δφ).

  This handles arbitrary readout placement — the readout drive does not
  need to sit on resonance with either parity branch.

- The resonator-basis loading approximation gives

      Δf_r / f_r ≈ - ΔC_Q / (2 C_r)                     (linear, low loading)
      ⇒  ΔC_Q [F]  =  -2 · C_r · Δf_r / f_r.

  `frequency_shift_to_quantum_capacitance` implements this directly.

- The loaded linewidth κ used by the phase helpers can be derived from
  internal and coupling quality factors via

      1/Q_loaded = 1/Q_i + 1/Q_c,        κ = f_r / Q_loaded.

  `kappa_from_quality_factors` returns κ given (f_r, Q_i, Q_c).
"""

import numpy as np


def phase_shift_to_frequency_shift(delta_phase_rad,
                                   kappa_hz,
                                   drive_detuning_hz=0.0):
    """
    Convert a measured cavity phase shift to a resonator frequency shift.

    Inverts the driven-Lorentzian phase response in closed form, so the
    readout drive does not need to be parked on the bare resonance.

    Parameters
    ----------
    delta_phase_rad : float or array_like
        Measured phase shift Δφ (radians) caused by the qubit-induced
        cavity frequency shift.
    kappa_hz : float
        Loaded resonator linewidth κ (FWHM) in Hz.
    drive_detuning_hz : float or array_like, optional
        Drive-vs-bare-resonator detuning Δ_d = f_drive - f_r [Hz].
        Default 0 (drive on resonance).

    Returns
    -------
    delta_f_hz : float or ndarray
        Implied resonator frequency shift Δf_r [Hz].
    """
    phi0 = np.arctan(2.0 * np.asarray(drive_detuning_hz) / kappa_hz)
    return (np.asarray(drive_detuning_hz)
            - 0.5 * kappa_hz * np.tan(phi0 + np.asarray(delta_phase_rad)))


def kappa_from_quality_factors(f_r_hz, q_i, q_c):
    """
    Loaded resonator linewidth κ from internal and coupling Q factors.

    Uses the parallel-Q rule

        1/Q_loaded = 1/Q_i + 1/Q_c,
        κ = f_r / Q_loaded
          = f_r · (1/Q_i + 1/Q_c).

    Convention: Q_c is taken to be real (symmetric hanger / notch
    geometry). For asymmetric hangers where the coupling Q is complex,
    pass `Re(1/Q_c)` × Q_c equivalent or use the loaded Q directly.

    Parameters
    ----------
    f_r_hz : float or array_like
        Resonator frequency f_r [Hz].
    q_i : float or array_like
        Internal (intrinsic) quality factor.
    q_c : float or array_like
        Coupling quality factor.

    Returns
    -------
    kappa_hz : float or ndarray
        Loaded resonator linewidth κ (FWHM) in Hz.
    """
    q_inv = 1.0 / np.asarray(q_i) + 1.0 / np.asarray(q_c)
    return np.asarray(f_r_hz) * q_inv


def frequency_shift_to_quantum_capacitance(delta_f_hz, f_r_hz, c_r_f):
    """
    Map a resonator frequency shift to a change in quantum capacitance.

    Uses the linear, low-loading approximation Δf_r/f_r ≈ -ΔC_Q/(2 C_r),
    valid when |ΔC_Q| ≪ C_r.

    Parameters
    ----------
    delta_f_hz : float or array_like
        Resonator frequency shift Δf_r [Hz]. Sign convention: positive
        Δf_r corresponds to negative ΔC_Q.
    f_r_hz : float
        Bare (unshifted) resonator frequency [Hz].
    c_r_f : float
        Effective resonator capacitance [F] that sets the loading scale.

    Returns
    -------
    delta_cq_f : float or ndarray
        Change in quantum capacitance ΔC_Q [F].
    """
    return -2.0 * c_r_f * np.asarray(delta_f_hz) / f_r_hz
