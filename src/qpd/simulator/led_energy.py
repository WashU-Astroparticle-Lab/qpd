"""Single-channel energy-to-tunneling model for a complete LED pulse."""

from __future__ import annotations

from dataclasses import dataclass, fields
import math

import numpy as np


def diffusion_trapping_probability(
    absorber_length_um: float, trap_length_um: float, *,
    diffusion_length_um: float = 300.0, scattering_time_us: float = 0.2,
    trap_diffusion_cm2_s: float = 2.0,
) -> float:
    """Appendix B, Eq. B1 of Ramanathan et al. (doi:10.1103/kqd2-spb1).

    absorber_length_um is the full width w in Fig. 11. The scattering to
    dwell time ratio uses the trap length l and dwell time l²/D.
    """
    for value in (absorber_length_um, trap_length_um, diffusion_length_um,
                  scattering_time_us, trap_diffusion_cm2_s):
        if not math.isfinite(value) or value <= 0:
            raise ValueError("diffusion parameters must be positive and finite")
    alpha = absorber_length_um / diffusion_length_um
    beta = scattering_time_us * trap_diffusion_cm2_s * 100 / trap_length_um**2
    # B1 simplifies to 1 / (cosh(alpha/2) + beta * sinh(alpha/2)).
    return 2 * math.exp(-alpha / 2) / (
        1 + math.exp(-alpha) + beta * -math.expm1(-alpha)
    )


@dataclass
class LEDPopulationTruth:
    """One entry per pulse, in rising-edge order. Mobile yield is effective."""

    absorber_qp: np.ndarray
    captured_qp: np.ndarray
    mobile_qp: np.ndarray
    lambda_ch: np.ndarray


@dataclass(frozen=True)
class LEDEnergyModel:
    """Fixed substrate deposit per pulse and one channel's effective collection.

    eta_ce includes the channel share and absorber pair breaking. Fixed
    multiplication_gain converts captures to an effective mobile population.
    The defaults describe an Al absorber and two equal Hf OCS traps.
    The mobile population is shared equally across both 500 µm³ volumes.
    """

    deposited_energy_ev: float
    eta_ce: float
    absorber_gap_ev: float = 190e-6
    fano_factor: float = 0.2
    trapping_probability: float = 0.2
    multiplication_gain: float = 2.375  # 0.5 * (190 µeV / 40 µeV)
    tunneling_constant_um3_s: float = 20e3
    qp_lifetime_s: float = 1e-3
    trap_volume_um3: float = 1000.0  # combined volume for the channel population

    def __post_init__(self) -> None:
        for field in fields(self):
            value = getattr(self, field.name)
            if not math.isfinite(value) or value < 0:
                raise ValueError(f"{field.name} must be finite and nonnegative")
        if self.absorber_gap_ev == 0 or self.trap_volume_um3 == 0:
            raise ValueError("absorber gap and trap volume must be positive")
        if max(self.eta_ce, self.fano_factor, self.trapping_probability) > 1:
            raise ValueError("efficiencies and absorber Fano factor must be in [0, 1]")

    def sample(self, n_pulses: int, rng: np.random.Generator) -> LEDPopulationTruth:
        """Draw integer absorber and captured counts, then the expected flips."""
        mean = self.deposited_energy_ev * self.eta_ce / self.absorber_gap_ev
        absorber = _sample_absorber_counts(mean, self.fano_factor, n_pulses, rng)
        captured = rng.binomial(absorber, self.trapping_probability)
        mobile = captured * self.multiplication_gain
        crossings_per_qp = (
            self.tunneling_constant_um3_s * self.qp_lifetime_s / self.trap_volume_um3
        )
        return LEDPopulationTruth(absorber, captured, mobile, crossings_per_qp * mobile)


def _sample_absorber_counts(mean, fano, size, rng):
    """Moment-matched integer approximation to the absorber production cascade.

    Start with floor(mean) trials of probability p and one fractional trial
    of probability frac(mean)*p. An independent Poisson contribution restores
    the mean. The variance is mean - p²*(floor(mean) + frac(mean)²).
    At p=1 this reaches the minimum possible variance for an integer count.
    """
    if mean == 0 or fano == 1:
        return rng.poisson(mean, size)
    whole = math.floor(mean)
    fraction = mean - whole
    if mean * fano <= fraction * (1 - fraction):
        return whole + rng.binomial(1, fraction, size)
    p = math.sqrt(mean * (1 - fano) / (whole + fraction**2))
    return (rng.binomial(whole, p, size) + rng.binomial(1, fraction * p, size)
            + rng.poisson(mean * (1 - p), size))
