"""Stationary dressed-state S21 forward model + lock-in bandwidth (sinc) window.

Demonstrates the two simulator upgrades at the WashU QPD device point:

  1. ``QPD.compute_s21_dressed`` -- complex S21(f) from a full-Hamiltonian
     diagonalization under the stationary-population approximation, in the Probst
     Eq. (1) notch convention (so it overlays on resonator_tools circle fits).
  2. ``apply_lockin_window`` -- the finite-bandwidth lock-in measurement, a sinc
     window of first null at ``df`` (boxcar integration over T = 1/df), convolved
     onto the complex S21.

Run: ``python examples/example_readout_window.py``  (writes example_readout_window.png).
"""
import numpy as np
import matplotlib.pyplot as plt

from qpd import QPD
from qpd.theory.driven import apply_lockin_window, measured_s21_dressed

# --- WashU QPD device point -------------------------------------------------
E_J, E_C = 6.95e9, 0.695e9
F_R, G = 6.936e9, 40e6
KAPPA, KAPPA_C = 115.70e3, 69.82e3      # total / coupling-port rates [Hz]


def main():
    qpd = QPD(e_j_hz=E_J, e_c_hz=E_C)

    # Probe window around the dressed cavity (chi ~ +0.47 MHz at this point).
    # Span is wide enough that the widest sinc window (df=200 kHz, 20 lobes) fits.
    freqs = np.linspace(F_R - 5e6, F_R + 7e6, 12001)

    # --- Parity-resolved true S21 (no measurement bandwidth) ----------------
    s21 = {
        p: qpd.compute_s21_dressed(
            0.0, G, F_R, freqs, parity=p, kappa_hz=KAPPA, kappa_c_hz=KAPPA_C,
        )
        for p in ("even", "odd")
    }

    # --- Lock-in bandwidth sweep on the even-parity trace -------------------
    df_list = [5e3, 50e3, 200e3]            # lock-in resolutions to compare
    windowed = {
        df: apply_lockin_window(freqs, s21["even"], df_hz=df) for df in df_list
    }

    # --- Report: true chi vs measured notch depth/width vs df ---------------
    def notch_depth(s):                      # peak deviation from baseline (=1)
        return float(np.max(np.abs(s - 1.0)))

    def fwhm_hz(s):                          # full width at half max of |S21-1|
        dev = np.abs(s - 1.0)
        half = 0.5 * dev.max()
        above = np.where(dev >= half)[0]
        return float(freqs[above[-1]] - freqs[above[0]]) if above.size else np.nan

    f_peak = freqs[np.argmax(np.abs(s21["even"] - 1.0))]
    print(f"WashU point: f_r = {F_R/1e9:.6f} GHz, kappa = {KAPPA/1e3:.1f} kHz")
    print(f"True (even) chi = f_peak - f_r = {(f_peak - F_R)/1e3:+.1f} kHz")
    print(f"True (even) notch depth = {notch_depth(s21['even']):.4f}, "
          f"FWHM = {fwhm_hz(s21['even'])/1e3:.1f} kHz")
    print("-" * 56)
    print(f"{'df [kHz]':>10} {'depth':>10} {'FWHM [kHz]':>12}")
    print(f"{'(true)':>10} {notch_depth(s21['even']):>10.4f} "
          f"{fwhm_hz(s21['even'])/1e3:>12.1f}")
    for df in df_list:
        print(f"{df/1e3:>10.0f} {notch_depth(windowed[df]):>10.4f} "
              f"{fwhm_hz(windowed[df])/1e3:>12.1f}")
    print("-" * 56)
    print("Note: df << kappa is faithful; df >~ kappa fills the notch and")
    print("broadens it (apparent Q drops). The symmetric window does not move")
    print("a symmetric peak -- a chi bias appears only for asymmetric lines.")

    # --- Figure -------------------------------------------------------------
    with plt.style.context(QPD._style_path):
        fig, ax = plt.subplots(2, 2, figsize=(9, 6), tight_layout=True)
        fghz = freqs / 1e9

        # (top-left) parity-resolved magnitude
        for p, c in (("even", "C0"), ("odd", "C3")):
            ax[0, 0].plot(fghz, np.abs(s21[p]), color=c, label=f"{p}")
        ax[0, 0].set(ylabel="|S21|", title="Parity-resolved true S21")
        ax[0, 0].legend()

        # (top-right) parity-resolved phase
        for p, c in (("even", "C0"), ("odd", "C3")):
            ax[0, 1].plot(fghz, np.angle(s21[p]), color=c, label=f"{p}")
        ax[0, 1].set(ylabel="arg S21 [rad]", title="Phase")

        # (bottom-left) bandwidth window on |S21| (even)
        ax[1, 0].plot(fghz, np.abs(s21["even"]), "k", lw=1.5, label="true")
        for df in df_list:
            ax[1, 0].plot(fghz, np.abs(windowed[df]), label=f"df={df/1e3:.0f} kHz")
        ax[1, 0].set(xlabel="frequency [GHz]", ylabel="|S21|",
                     title="Lock-in sinc window (even)")
        ax[1, 0].legend()

        # (bottom-right) complex-plane resonance circle
        ax[1, 1].plot(s21["even"].real, s21["even"].imag, "k", lw=1.5, label="true")
        for df in df_list:
            ax[1, 1].plot(windowed[df].real, windowed[df].imag,
                          label=f"df={df/1e3:.0f} kHz")
        ax[1, 1].set(xlabel="Re S21", ylabel="Im S21",
                     title="Resonance circle")
        ax[1, 1].set_aspect("equal", "box")
        ax[1, 1].legend()

        out = "example_readout_window.png"
        fig.savefig(out, dpi=150)
        print(f"\nSaved figure -> {out}")
        plt.show()


if __name__ == "__main__":
    main()
