"""Quasiparticle burst model: EMG-profiled bursts of QP tunneling events.

A *burst* is an energy-deposition event that injects a number of quasiparticle
(QP) tunneling events whose arrival times, relative to the burst onset, follow
an exponentially modified Gaussian (EMG) distribution: the sum of a
``Normal(mu, sigma)`` and an independent ``Exponential(tau)`` random variable.

Each tunneling event flips the transmon charge parity, exactly like a
background telegraph switch (see :mod:`qpd.simulator.parity`). The bursts are
treated as an independent point process superimposed on the background CTMC;
``VNASimulator`` merges the two flip-time sets when building the parity
trajectory. All times are in SECONDS.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class BurstTruth:
    """Ground-truth record for a single quasiparticle burst.

    Attributes
    ----------
    t_arrival : float
        Burst onset time (the EMG location reference ``t_burst``), seconds.
    n_qp : int
        Actual number of QP tunneling events drawn for this burst.
    t_start : float
        Time of the earliest event in the burst (min of ``event_times``),
        seconds. Equals ``t_arrival`` when ``n_qp == 0`` (degenerate span).
    t_end : float
        Time of the latest event in the burst (max of ``event_times``),
        seconds. Equals ``t_arrival`` when ``n_qp == 0``.
    event_times : np.ndarray
        Sorted absolute event times for this burst, shape ``(n_qp,)``, seconds.
    """

    t_arrival: float
    n_qp: int
    t_start: float
    t_end: float
    event_times: np.ndarray


@dataclass
class QuasiparticleBurstModel:
    """Generate EMG-profiled bursts of quasiparticle tunneling events.

    Parameters
    ----------
    times : np.ndarray
        Burst arrival (onset) times in seconds, shape ``(n_bursts,)``. Use an
        explicit list for a controlled demo, or :func:`poisson_burst_times` to
        draw random arrivals.
    tau : float
        Exponential time constant of the EMG, seconds (e.g. ``3.7e-3``).
    mu : float
        Gaussian mean offset of the EMG, seconds (e.g. ``1.2e-3``).
    sigma : float
        Gaussian standard deviation of the EMG, seconds (e.g. ``0.4e-3``).
    expected_n_qp : float | np.ndarray, default 10.0
        Expected number of QP tunneling events per burst (Poisson mean).
        Scalar (shared by all bursts) or an array broadcastable to ``times``.
        The actual per-burst count is drawn ``N ~ Poisson(expected_n_qp)`` and
        recorded as truth.

    Notes
    -----
    Per burst ``i``: draw ``N_i ~ Poisson(expected_n_qp_i)``, then ``N_i``
    offsets ``= rng.normal(mu, sigma, N_i) + rng.exponential(tau, N_i)``; the
    absolute event times are ``times[i] + offsets``. Offsets can be negative
    (the Gaussian lower tail), so a burst may emit a few events slightly before
    its onset; ``BurstTruth.t_start`` reflects the actual earliest event.
    """

    times: np.ndarray
    tau: float
    mu: float
    sigma: float
    expected_n_qp: float | np.ndarray = 10.0

    def __post_init__(self) -> None:
        self.times = np.asarray(self.times, dtype=float)
        if self.times.ndim != 1:
            raise ValueError("times must be 1-D")
        if self.tau < 0 or self.sigma < 0:
            raise ValueError("tau and sigma must be non-negative")
        exp = np.asarray(self.expected_n_qp, dtype=float)
        if exp.ndim == 0:
            exp = np.full(self.times.shape, float(exp))
        elif exp.shape != self.times.shape:
            raise ValueError(
                "expected_n_qp must be scalar or match the shape of times"
            )
        if np.any(exp < 0):
            raise ValueError("expected_n_qp must be non-negative")
        self._expected = exp

    def sample(
        self, rng: np.random.Generator
    ) -> tuple[np.ndarray, list[BurstTruth]]:
        """Sample all burst event times and per-burst truth.

        Parameters
        ----------
        rng : np.random.Generator
            Random generator (seeded upstream by ``VNASimulator.simulate``).

        Returns
        -------
        event_times : np.ndarray
            Sorted absolute flip times across ALL bursts, shape
            ``(sum_i N_i,)``, seconds. An empty float array if no events are
            produced.
        truth : list[BurstTruth]
            One :class:`BurstTruth` per burst, in the same order as ``times``.
        """
        truth: list[BurstTruth] = []
        all_events: list[np.ndarray] = []
        for t0, lam in zip(self.times, self._expected):
            t0 = float(t0)
            n = int(rng.poisson(lam))
            if n == 0:
                truth.append(
                    BurstTruth(
                        t_arrival=t0,
                        n_qp=0,
                        t_start=t0,
                        t_end=t0,
                        event_times=np.empty(0, dtype=float),
                    )
                )
                continue
            offsets = rng.normal(self.mu, self.sigma, n) + rng.exponential(
                self.tau, n
            )
            ev = np.sort(t0 + offsets)
            all_events.append(ev)
            truth.append(
                BurstTruth(
                    t_arrival=t0,
                    n_qp=n,
                    t_start=float(ev[0]),
                    t_end=float(ev[-1]),
                    event_times=ev,
                )
            )
        if all_events:
            event_times = np.sort(np.concatenate(all_events))
        else:
            event_times = np.empty(0, dtype=float)
        return event_times, truth


def poisson_burst_times(
    rate_hz: float,
    duration: float,
    rng: np.random.Generator,
    t_start: float = 0.0,
) -> np.ndarray:
    """Sample Poisson-random burst arrival times over ``[t_start, duration]``.

    A helper for the "random arrivals" mode: feed the result as the ``times``
    of a :class:`QuasiparticleBurstModel` so explicit and random modes share
    one model.

    Parameters
    ----------
    rate_hz : float
        Mean burst arrival rate, Hz. The expected number of bursts is
        ``rate_hz * (duration - t_start)``.
    duration : float
        End of the arrival window, seconds.
    rng : np.random.Generator
        Random generator.
    t_start : float, default 0.0
        Start of the arrival window, seconds.

    Returns
    -------
    np.ndarray
        Sorted burst arrival times in ``[t_start, duration]``, seconds.
    """
    if rate_hz < 0:
        raise ValueError("rate_hz must be non-negative")
    window = float(duration) - float(t_start)
    if window <= 0 or rate_hz == 0:
        return np.empty(0, dtype=float)
    n = int(rng.poisson(rate_hz * window))
    if n == 0:
        return np.empty(0, dtype=float)
    return np.sort(rng.uniform(t_start, duration, n))
