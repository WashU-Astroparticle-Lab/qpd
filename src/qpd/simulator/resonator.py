"""Notch-type resonator S21 model (Probst et al., RSI 86, 024706 (2015), Eq. 1)."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Union

import numpy as np

ArrayLike = Union[float, np.ndarray]


@dataclass
class ResonatorConfig:
    """Notch-geometry resonator + measurement-environment parameters."""

    f_r: float
    q_i: float
    q_c_abs: float
    phi: float = 0.0
    a: float = 1.0
    alpha: float = 0.0
    tau: float = 0.0

    def q_c_complex(self) -> complex:
        return self.q_c_abs * np.exp(-1j * self.phi)

    def q_l(self) -> float:
        inv_q_l = 1.0 / self.q_i + np.real(1.0 / self.q_c_complex())
        return 1.0 / inv_q_l


def notch_s21(
    f_drive: ArrayLike,
    f_r: ArrayLike,
    cfg: ResonatorConfig,
) -> np.ndarray:
    """Probst Eq. 1.

    S21(f) = a e^{iα} e^{-2π i f τ}
             · [ 1 - (Q_l / |Q_c|) e^{iφ} / (1 + 2 i Q_l (f / f_r - 1)) ]

    Either `f_drive` or `f_r` (or both) can be arrays; the result broadcasts.
    """
    f_drive = np.asarray(f_drive, dtype=float)
    f_r_arr = np.asarray(f_r, dtype=float)

    inv_q_l = 1.0 / cfg.q_i + np.real(1.0 / cfg.q_c_complex())
    q_l = 1.0 / inv_q_l

    ideal = 1.0 - (q_l / cfg.q_c_abs) * np.exp(1j * cfg.phi) / (
        1.0 + 2j * q_l * (f_drive / f_r_arr - 1.0)
    )
    env = cfg.a * np.exp(1j * cfg.alpha) * np.exp(-2j * np.pi * f_drive * cfg.tau)
    return env * ideal
