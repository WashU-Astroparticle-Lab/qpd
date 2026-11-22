"""Fixed-device calibration and shared simulator context for LED datasets."""
from __future__ import annotations

from dataclasses import asdict, dataclass, field

import numpy as np
from scipy.optimize import minimize_scalar

from qpd import (
    LEDPulseModel, QPD, ResonatorConfig, VNASimulator, WhiteGaussianNoise,
    led_count_sweep, notch_s21,
)


SPLITS = ("train", "validation", "test")
@dataclass(frozen=True)
class LEDDatasetConfig:
    """Nominal public settings; random seeds are deliberately not configuration."""

    snr_max_values: tuple[float, ...] = field(
        default_factory=lambda: tuple(np.geomspace(0.5, 5.0, 13)))
    mean_tunnels_values: tuple[float, ...] = field(
        default_factory=lambda: tuple(led_count_sweep()))
    chunks_per_split: dict[str, int] = field(
        default_factory=lambda: {"train": 4, "validation": 2, "test": 2})
    duration_s: float = 5.0
    sample_rate_hz: float = 100_000.0
    e_j_hz: float = 6.5e9
    e_c_hz: float = 0.321e9
    coupling_g_hz: float = 150e6
    f_r_hz: float = 6.936e9
    # Even-branch resonance at n_g=0.25 plus 1 kHz for the default device.
    f_drive_hz: float = 6.939533642344854e9
    q_i: float = 1e5
    q_c_abs: float = 1.2e5
    phi: float = 0.03
    amplitude: float = 0.5
    alpha: float = 0.4
    cable_delay_s: float = 50e-9
    on_time_s: float = 50e-6
    off_time_s: float = 0.300
    first_pulse_time_s: float = 0.050
    burst_mu_s: float = 1.2e-3
    burst_sigma_s: float = 0.4e-3
    burst_tau_s: float = 3.7e-3
    quiescent_rate_hz: float = 20.0
    charge_jump_rate_hz: float = 1 / 60
    n_g_grid_points: int = 401
    charge_cutoff: int = 30

    def __post_init__(self):
        for name, positive in (("snr_max_values", True), ("mean_tunnels_values", False)):
            values = np.asarray(getattr(self, name), dtype=float)
            if (values.ndim != 1 or not values.size or not np.all(np.isfinite(values))
                    or np.any(values <= 0 if positive else values < 0)
                    or np.any(np.diff(values) <= 0)):
                raise ValueError(f"{name} must be a nonempty, increasing finite grid")
            object.__setattr__(self, name, tuple(float(v) for v in values))
        if set(self.chunks_per_split) != set(SPLITS):
            raise ValueError(f"chunks_per_split must contain exactly {SPLITS}")
        for split, count in self.chunks_per_split.items():
            if isinstance(count, bool) or not isinstance(count, (int, np.integer)) or count < 1:
                raise ValueError(f"{split} chunks per condition must be a positive integer")
        object.__setattr__(self, "chunks_per_split", {
            key: int(self.chunks_per_split[key]) for key in SPLITS})
        positive = ("duration_s", "sample_rate_hz", "e_j_hz", "e_c_hz", "coupling_g_hz",
                    "f_r_hz", "f_drive_hz", "q_i", "q_c_abs", "amplitude", "on_time_s")
        nonnegative = ("off_time_s", "first_pulse_time_s", "burst_mu_s", "burst_sigma_s",
                       "burst_tau_s", "quiescent_rate_hz", "charge_jump_rate_hz")
        for name in positive + nonnegative + ("phi", "alpha", "cable_delay_s"):
            value = getattr(self, name)
            if (np.ndim(value) != 0 or not np.isfinite(value)
                    or (name in positive and value <= 0)
                    or (name in nonnegative and value < 0)):
                raise ValueError(f"invalid {name}")
        for name, minimum in (("n_g_grid_points", 3), ("charge_cutoff", 2)):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, (int, np.integer)) or value < minimum:
                raise ValueError(f"{name} must be an integer >= {minimum}")
        if self.n_samples < 2:
            raise ValueError("duration_s * sample_rate_hz must produce at least two samples")

    @property
    def n_samples(self) -> int:
        return int(round(self.duration_s * self.sample_rate_hz))

    def to_dict(self) -> dict:
        return asdict(self)


class _CachedLEDSimulator(VNASimulator):
    """Use one full-charge interpolation table for every dataset chunk."""

    def _chi_grid(self, n_g_lo, n_g_hi):
        return self.cached_chi


@dataclass
class LEDSimulationContext:
    config: LEDDatasetConfig
    qpd: QPD
    resonator: ResonatorConfig
    chi_grid: tuple[np.ndarray, np.ndarray, np.ndarray]
    max_separation: float

    def branch_response(self, n_g):
        grid, even, odd = self.chi_grid
        return tuple(notch_s21(
            self.config.f_drive_hz,
            self.resonator.f_r + np.interp(n_g, grid, values), self.resonator,
        ) for values in (even, odd))

    def noise_sigma(self, snr_max: float) -> float:
        if np.ndim(snr_max) != 0 or not np.isfinite(snr_max) or snr_max <= 0:
            raise ValueError("snr_max must be finite and positive")
        return self.max_separation / float(snr_max)


def build_led_context(config: LEDDatasetConfig | None = None) -> LEDSimulationContext:
    """Diagonalize once; calibrate the maximum of the same interpolated model."""
    config = config or LEDDatasetConfig()
    qpd = QPD(e_j_hz=config.e_j_hz, e_c_hz=config.e_c_hz,
              coupling_g_hz=config.coupling_g_hz)
    resonator = ResonatorConfig(
        f_r=config.f_r_hz, q_i=config.q_i, q_c_abs=config.q_c_abs,
        phi=config.phi, a=config.amplitude, alpha=config.alpha, tau=config.cable_delay_s)
    simulator = VNASimulator(
        qpd=qpd, resonator=resonator, f_drive=config.f_drive_hz,
        sample_rate=config.sample_rate_hz, n_g_grid_points=config.n_g_grid_points,
        charge_cutoff=config.charge_cutoff)
    chi_grid = simulator._chi_grid(-0.5, 0.5)
    context = LEDSimulationContext(config, qpd, resonator, chi_grid, 0.0)

    def separation(charge):
        even, odd = context.branch_response(charge)
        return np.abs(even - odd)

    # Locate all maxima densely, then refine each on the actual interpolated
    # readout curve. The physical response is periodic at the half-open edges.
    charge = np.linspace(-0.5, np.nextafter(0.5, -np.inf), 16385)
    distances = separation(charge)
    if not np.all(np.isfinite(distances)):
        raise ValueError("non-finite response in SNR calibration")
    maxima = np.flatnonzero((distances[1:-1] >= distances[:-2])
                           & (distances[1:-1] >= distances[2:])) + 1
    maximum = float(np.max(distances))
    for index in maxima:
        optimum = minimize_scalar(
            lambda x: -float(separation(x)),
            bounds=(charge[index - 1], charge[index + 1]), method="bounded",
            options={"xatol": 1e-13})
        maximum = max(maximum, -float(optimum.fun))
    if not np.isfinite(maximum) or maximum <= 0:
        raise ValueError("SNR calibration requires nonzero finite parity separation")
    context.max_separation = maximum
    for values in chi_grid:
        values.flags.writeable = False
    return context
