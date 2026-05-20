"""Pluggable noise models for VNA I/Q traces."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol, runtime_checkable

import numpy as np


@runtime_checkable
class NoiseModel(Protocol):
    """Additive complex (I + iQ) noise sampler.

    Implementations return a length-`n_samples` complex array to be added
    to the noiseless S21 trace. `dt` is the sample spacing in seconds —
    white-noise implementations ignore it, but spectrally-shaped models
    (e.g. a future ``PSDNoise``) use it to set the frequency grid.
    """

    def sample(
        self,
        n_samples: int,
        dt: float,
        rng: np.random.Generator,
    ) -> np.ndarray: ...


@dataclass
class WhiteGaussianNoise:
    """Independent Gaussian noise on I and Q with per-quadrature std `sigma`."""

    sigma: float

    def sample(
        self,
        n_samples: int,
        dt: float,
        rng: np.random.Generator,
    ) -> np.ndarray:
        del dt
        return rng.normal(0.0, self.sigma, n_samples) + 1j * rng.normal(
            0.0, self.sigma, n_samples
        )
