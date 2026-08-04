"""Shared scenario builder for parity-flip reconstruction checks and studies.

The reference operating point is the one in ``notebooks/simulation.ipynb`` (the
QPD4 VNA fit, E_J/E_C = 12, g = 150 MHz, sigma = 2e-4, 100 kHz sampling), with
two changes that reflect how the measurement is actually run:

* the offset-charge sawtooth ramps at **500 Hz**, not 9 Hz. The ramp spans 5.5
  Cooper pairs = 11 half-pairs, so it resets by half a pair modulo one every
  2 ms -- and each reset swaps the two parity branches, which is
  indistinguishable from a parity flip event by event;
* the background tunnelling rate is **10 Hz** per state, so the ramp resets
  outnumber real flips by roughly 50:1.

Parity is a pure telegraph: it holds until the next quasiparticle tunnels and
each tunnel toggles it, with equal dwell times for the two states.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from qpd import (QPD, ChargeJumpEvents, ConstantNg, QuasiparticleBurstModel,
                 ResonatorConfig, SawtoothNg, VNASimulator,
                 WhiteGaussianNoise)

# QPD4 notch fit: dressed resonance and quality factors.
FR = 6489353036.703264
QI = 39669.674
QC = 81560.897
NG_MIN, NG_MAX = -0.5, 5.0
SPAN = NG_MAX - NG_MIN  # 5.5 Cooper pairs = 11 half-pairs (odd -> resets swap)
SAMPLE_RATE = 1e5
# The notebook parks the readout LO 15 MHz off the dressed line, which is what
# sets the per-sample SNR near unity. Kept so results are comparable.
BARE_SOLVE_OFFSET = 15e6


@dataclass
class Scenario:
    """A built simulator plus the truth quantities a reconstruction is scored on."""

    sim: VNASimulator
    duration: float
    fold_period: float  # [s] expected |chi| fold period = 0.5 / slope
    ramp_period: float  # [s] expected sawtooth reset period
    half_pairs_per_ramp: float  # ramp span in half Cooper pairs
    jump_times: np.ndarray = field(default_factory=lambda: np.empty(0))
    jump_deltas: np.ndarray = field(default_factory=lambda: np.empty(0))
    label: str = ""

    def simulate(self, seed: int = 0):
        return self.sim.simulate(duration=self.duration, seed=seed)


def build_scenario(
    ramp_hz: float = 500.0,
    tunnel_rate_hz: float = 10.0,
    duration: float = 5.0,
    sigma: float = 2e-4,
    e_j_hz: float = 8.335e9,
    e_c_hz: float = 0.695e9,
    g_hz: float = 150e6,
    burst_rate_hz: float = 1.0,
    burst_n_qp: float | np.ndarray = 15.0,
    burst_times: np.ndarray | None = None,
    charge_jump_rate_hz: float = 0.0,
    charge_jumps: tuple[np.ndarray, np.ndarray] | None = None,
    sample_rate: float = SAMPLE_RATE,
    seed: int = 0,
    ramp_t0: float = 0.0,
    label: str = "",
) -> Scenario:
    """Build the reference scenario with the knobs a study wants to vary.

    ``charge_jumps`` pins the jump times and amplitudes explicitly (used by the
    regression gate, where the worst case matters); otherwise jumps are drawn at
    ``charge_jump_rate_hz`` with amplitudes uniform on (-0.5, 0.5).
    """
    qpd = QPD(e_j_hz=e_j_hz, e_c_hz=e_c_hz)
    qpd.coupling_g_hz = g_hz
    # Back out the bare cavity so the qubit pull places the dressed even line
    # at the measured FR, offset as in the notebook.
    f_bare = FR
    for _ in range(6):
        _, chi = qpd.compute_dispersive_matrix(0.0, g_hz, f_bare,
                                               num_levels=2, parity="even")
        f_bare = FR - chi[0] + BARE_SOLVE_OFFSET
    resonator = ResonatorConfig(f_r=f_bare, q_i=QI, q_c_abs=QC,
                                phi=0.03, a=0.5, alpha=0.4, tau=50e-9)

    rng = np.random.default_rng(seed + 991)
    slope = SPAN * ramp_hz
    offset_charge = SawtoothNg(NG_MIN, NG_MAX, slope=slope, t0=ramp_t0)

    if charge_jumps is not None:
        jump_t = np.asarray(charge_jumps[0], dtype=float)
        jump_d = np.asarray(charge_jumps[1], dtype=float)
    else:
        n_jump = rng.poisson(charge_jump_rate_hz * duration)
        jump_t = np.sort(rng.uniform(0.0, duration, n_jump))
        jump_d = rng.uniform(-0.5, 0.5, n_jump)
    if jump_t.size:
        offset_charge = offset_charge + ChargeJumpEvents(times=jump_t,
                                                         deltas=jump_d)

    bursts = None
    if burst_times is not None:
        onsets = np.asarray(burst_times, dtype=float)
    elif burst_rate_hz > 0:
        onsets = np.sort(rng.uniform(0.0, duration,
                                     rng.poisson(burst_rate_hz * duration)))
    else:
        onsets = np.empty(0)
    if onsets.size:
        n_qp = burst_n_qp
        if np.ndim(n_qp) > 0:
            n_qp = np.resize(np.asarray(n_qp, dtype=float), onsets.size)
        bursts = QuasiparticleBurstModel(times=onsets, tau=3.7e-3, mu=1.2e-3,
                                         sigma=0.4e-3, expected_n_qp=n_qp)

    sim = VNASimulator(
        qpd=qpd, resonator=resonator, f_drive=FR + 1e3,
        sample_rate=sample_rate,
        gamma_even_to_odd=tunnel_rate_hz, gamma_odd_to_even=tunnel_rate_hz,
        noise=WhiteGaussianNoise(sigma=sigma),
        offset_charge=offset_charge, quasiparticle_bursts=bursts,
    )
    return Scenario(sim=sim, duration=duration, fold_period=0.5 / slope,
                    ramp_period=1.0 / ramp_hz, half_pairs_per_ramp=SPAN / 0.5,
                    jump_times=jump_t, jump_deltas=jump_d, label=label)


def build_static_scenario(
    n_g: float = 0.0,
    tunnel_rate_hz: float = 10.0,
    duration: float = 5.0,
    sigma: float = 2e-4,
    e_j_hz: float = 8.335e9,
    e_c_hz: float = 0.695e9,
    g_hz: float = 150e6,
    burst_rate_hz: float = 0.0,
    burst_n_qp: float | np.ndarray = 15.0,
    sample_rate: float = SAMPLE_RATE,
    seed: int = 0,
    label: str = "",
) -> Scenario:
    """Constant-offset-charge scenario: two stationary blobs, telegraph hopping.

    This is the regime :func:`qpd.reconstruction.reconstruct_parity_flips_static`
    targets. There is no ramp, so no fold period, no reset comb and no
    parity-blind point sweeping past -- but the contrast is fixed by ``n_g``.
    Best contrast is at ``n_g = 0`` (mod 0.5); ``n_g = 0.25`` (mod 0.5) is the
    parity-blind charge, where the two branches coincide and no algorithm can
    recover the flips.
    """
    qpd = QPD(e_j_hz=e_j_hz, e_c_hz=e_c_hz)
    qpd.coupling_g_hz = g_hz
    f_bare = FR
    for _ in range(6):
        _, chi = qpd.compute_dispersive_matrix(0.0, g_hz, f_bare,
                                               num_levels=2, parity="even")
        f_bare = FR - chi[0] + BARE_SOLVE_OFFSET
    resonator = ResonatorConfig(f_r=f_bare, q_i=QI, q_c_abs=QC,
                                phi=0.03, a=0.5, alpha=0.4, tau=50e-9)

    rng = np.random.default_rng(seed + 991)
    bursts = None
    if burst_rate_hz > 0:
        onsets = np.sort(rng.uniform(0.0, duration,
                                     rng.poisson(burst_rate_hz * duration)))
        if onsets.size:
            bursts = QuasiparticleBurstModel(times=onsets, tau=3.7e-3,
                                             mu=1.2e-3, sigma=0.4e-3,
                                             expected_n_qp=burst_n_qp)

    sim = VNASimulator(
        qpd=qpd, resonator=resonator, f_drive=FR + 1e3,
        sample_rate=sample_rate,
        gamma_even_to_odd=tunnel_rate_hz, gamma_odd_to_even=tunnel_rate_hz,
        noise=WhiteGaussianNoise(sigma=sigma),
        offset_charge=ConstantNg(n_g), quasiparticle_bursts=bursts,
    )
    return Scenario(sim=sim, duration=duration, fold_period=float("nan"),
                    ramp_period=float("nan"), half_pairs_per_ramp=float("nan"),
                    label=label or f"static n_g={n_g}")
