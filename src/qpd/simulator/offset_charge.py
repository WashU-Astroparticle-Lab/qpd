"""Composable offset-charge n_g(t) models."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Sequence

import numpy as np


class OffsetChargeModel:
    """Base class. Subclasses implement ``evaluate(t)``.

    Any two models can be summed via ``a + b`` to produce a `CompositeNg`.
    This lets a sawtooth drift be combined with discrete charge-jump events,
    or any other future component, without ceremony.
    """

    def evaluate(self, t: np.ndarray) -> np.ndarray:
        raise NotImplementedError

    def __add__(self, other: "OffsetChargeModel") -> "CompositeNg":
        if not isinstance(other, OffsetChargeModel):
            return NotImplemented
        return CompositeNg(self, other)


@dataclass
class ConstantNg(OffsetChargeModel):
    value: float

    def evaluate(self, t: np.ndarray) -> np.ndarray:
        return np.full_like(np.asarray(t, dtype=float), self.value)


@dataclass
class SawtoothNg(OffsetChargeModel):
    """Linear ramp at `slope` [1/s] between n_g_min and n_g_max.

    For positive slope, n_g(t0) = n_g_min and ramps upward, snapping back
    to n_g_min when it reaches n_g_max. For negative slope, n_g(t0) = n_g_max
    and ramps downward. This reproduces the slow charge-drift sawtooth in
    Fig. 7 (top) of arXiv:2405.17192.
    """

    n_g_min: float
    n_g_max: float
    slope: float
    t0: float = 0.0

    def evaluate(self, t: np.ndarray) -> np.ndarray:
        t = np.asarray(t, dtype=float)
        span = self.n_g_max - self.n_g_min
        if span <= 0:
            raise ValueError("SawtoothNg requires n_g_max > n_g_min")
        phase_raw = self.slope * (t - self.t0)
        if self.slope >= 0:
            return self.n_g_min + np.mod(phase_raw, span)
        return self.n_g_max + np.mod(phase_raw, -span)


@dataclass
class ChargeJumpEvents(OffsetChargeModel):
    """Step function: cumulative additive jumps at the listed times."""

    times: np.ndarray
    deltas: np.ndarray

    def __post_init__(self) -> None:
        self.times = np.asarray(self.times, dtype=float)
        self.deltas = np.asarray(self.deltas, dtype=float)
        if self.times.shape != self.deltas.shape:
            raise ValueError("times and deltas must have the same shape")
        if not np.all(np.diff(self.times) >= 0):
            raise ValueError("times must be monotonically non-decreasing")

    def evaluate(self, t: np.ndarray) -> np.ndarray:
        t = np.asarray(t, dtype=float)
        cum = np.concatenate(([0.0], np.cumsum(self.deltas)))
        idx = np.searchsorted(self.times, t, side="right")
        return cum[idx]


@dataclass
class CompositeNg(OffsetChargeModel):
    """Sum of any number of OffsetChargeModel components."""

    components: Sequence[OffsetChargeModel] = field(default_factory=tuple)

    def __init__(self, *components: OffsetChargeModel) -> None:
        flat: list[OffsetChargeModel] = []
        for c in components:
            if isinstance(c, CompositeNg):
                flat.extend(c.components)
            else:
                flat.append(c)
        self.components = tuple(flat)

    def evaluate(self, t: np.ndarray) -> np.ndarray:
        t = np.asarray(t, dtype=float)
        out = np.zeros_like(t)
        for c in self.components:
            out = out + c.evaluate(t)
        return out
