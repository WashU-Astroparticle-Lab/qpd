"""Top-level VNA S21 timeseries simulator."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from qpd.theory.transmon import QPD

from .noise import NoiseModel, WhiteGaussianNoise
from .offset_charge import ConstantNg, OffsetChargeModel
from .parity import EVEN, generate_parity_trajectory
from .resonator import ResonatorConfig, notch_s21


@dataclass
class SimResult:
    t: np.ndarray
    i: np.ndarray
    q: np.ndarray
    parity: np.ndarray
    n_g: np.ndarray

    @property
    def iq(self) -> np.ndarray:
        return self.i + 1j * self.q


def _default_noise() -> NoiseModel:
    return WhiteGaussianNoise(sigma=0.0)


def _default_ng() -> OffsetChargeModel:
    return ConstantNg(0.25)


@dataclass
class VNASimulator:
    """Generate I/Q timeseries of a QPD-coupled notch resonator.

    For each sample on a uniform time grid, the parity-resolved
    ground-state dispersive shift χ at the current n_g is added to the
    bare resonance frequency, the S21 model is evaluated at the fixed
    drive frequency, and additive noise is applied. χ(n_g) is precomputed
    on a fine grid and interpolated to keep long traces fast.
    """

    qpd: QPD
    resonator: ResonatorConfig
    f_drive: float
    sample_rate: float
    gamma_even_to_odd: float
    gamma_odd_to_even: float
    noise: NoiseModel = field(default_factory=_default_noise)
    offset_charge: OffsetChargeModel = field(default_factory=_default_ng)
    n_g_grid_points: int = 401
    charge_cutoff: int = 30

    def _chi_grid(self, n_g_lo: float, n_g_hi: float):
        # Margin to avoid edge-extrapolation when the trajectory grazes a bound.
        if n_g_hi - n_g_lo < 1e-9:
            n_g_lo -= 0.05
            n_g_hi += 0.05
        else:
            margin = 0.02 * (n_g_hi - n_g_lo)
            n_g_lo -= margin
            n_g_hi += margin

        grid = np.linspace(n_g_lo, n_g_hi, self.n_g_grid_points)
        chi_even = np.empty_like(grid)
        chi_odd = np.empty_like(grid)
        g = self.qpd.coupling_g_hz
        for k, n_g in enumerate(grid):
            _, chi_e = self.qpd.compute_dispersive_matrix(
                float(n_g),
                coupling_g_hz=g,
                readout_freq_hz=self.f_drive,
                num_levels=2,
                parity="even",
                charge_cutoff=self.charge_cutoff,
            )
            _, chi_o = self.qpd.compute_dispersive_matrix(
                float(n_g),
                coupling_g_hz=g,
                readout_freq_hz=self.f_drive,
                num_levels=2,
                parity="odd",
                charge_cutoff=self.charge_cutoff,
            )
            chi_even[k] = chi_e[0]
            chi_odd[k] = chi_o[0]
        return grid, chi_even, chi_odd

    def simulate(self, duration: float, seed: int | None = None) -> SimResult:
        rng = np.random.default_rng(seed)
        n = int(round(duration * self.sample_rate))
        if n < 2:
            raise ValueError("duration * sample_rate must produce >= 2 samples")
        dt = 1.0 / self.sample_rate
        t = np.arange(n) * dt

        n_g_t = self.offset_charge.evaluate(t)

        parity_t = generate_parity_trajectory(
            t,
            self.gamma_even_to_odd,
            self.gamma_odd_to_even,
            rng,
        )

        n_g_grid, chi_even, chi_odd = self._chi_grid(
            float(np.min(n_g_t)), float(np.max(n_g_t))
        )
        chi_even_t = np.interp(n_g_t, n_g_grid, chi_even)
        chi_odd_t = np.interp(n_g_t, n_g_grid, chi_odd)
        chi_t = np.where(parity_t == EVEN, chi_even_t, chi_odd_t)

        f_r_eff = self.resonator.f_r + chi_t
        s21 = notch_s21(self.f_drive, f_r_eff, self.resonator)
        s21 = s21 + self.noise.sample(n, dt, rng)

        return SimResult(
            t=t,
            i=s21.real,
            q=s21.imag,
            parity=parity_t,
            n_g=n_g_t,
        )
