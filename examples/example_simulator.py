"""End-to-end demo of qpd.simulator: VNA S21 timeseries of a QPD transmon."""

import numpy as np
import matplotlib.pyplot as plt

from qpd import (
    QPD,
    VNASimulator,
    ResonatorConfig,
    WhiteGaussianNoise,
    SawtoothNg,
    ChargeJumpEvents,
)


def main():
    qpd = QPD(e_j_hz=8.335e9, e_c_hz=0.695e9)
    qpd.coupling_g_hz = 150e6

    resonator = ResonatorConfig(
        f_r=6.0e9,
        q_i=1e5,
        q_c_abs=1e4,
        phi=0.03,
        a=0.5,
        alpha=0.4,
        tau=50e-9,
    )

    n_g_model = SawtoothNg(
        n_g_min=-0.5, n_g_max=0.5, slope=0.05
    ) + ChargeJumpEvents(
        times=np.array([1.2, 3.5, 7.8]),
        deltas=np.array([+0.3, -0.2, +0.4]),
    )

    gamma_e_to_o = 200.0
    gamma_o_to_e = 500.0
    p_odd_expected = gamma_e_to_o / (gamma_e_to_o + gamma_o_to_e)

    sim = VNASimulator(
        qpd=qpd,
        resonator=resonator,
        f_drive=6.0e9,
        sample_rate=1e6,
        gamma_even_to_odd=gamma_e_to_o,
        gamma_odd_to_even=gamma_o_to_e,
        noise=WhiteGaussianNoise(sigma=2e-4),
        offset_charge=n_g_model,
    )

    result = sim.simulate(duration=10.0, seed=0)

    frac_odd = float(np.mean(result.parity == 1))
    print(f"Empirical P(odd) = {frac_odd:.4f}")
    print(f"Predicted P(odd) = {p_odd_expected:.4f}")

    with plt.style.context(QPD._style_path):
        fig, axes = plt.subplots(3, 1, figsize=(7, 8), constrained_layout=True)

        zoom = result.t < 0.05
        axes[0].plot(result.t[zoom] * 1e3, result.i[zoom], label="I", lw=0.5)
        axes[0].plot(result.t[zoom] * 1e3, result.q[zoom], label="Q", lw=0.5)
        axes[0].set_xlabel("Time [ms]")
        axes[0].set_ylabel("Voltage [a.u.]")
        axes[0].set_title("I(t), Q(t) — first 50 ms")
        axes[0].legend(loc="best")

        sc = axes[1].scatter(
            result.i, result.q, c=result.parity, cmap="coolwarm",
            s=1, alpha=0.3,
        )
        axes[1].set_xlabel("I")
        axes[1].set_ylabel("Q")
        axes[1].set_aspect("equal", adjustable="datalim")
        axes[1].set_title("IQ scatter (color = parity)")
        plt.colorbar(sc, ax=axes[1], label="parity (0=even, 1=odd)")

        axes[2].plot(result.t, result.n_g, lw=0.5)
        axes[2].set_xlabel("Time [s]")
        axes[2].set_ylabel(r"$n_g$")
        axes[2].set_title(r"$n_g(t)$: sawtooth + discrete jumps")

        fig.savefig("example_simulator.png", dpi=150)
        print("Saved example_simulator.png")


if __name__ == "__main__":
    main()
