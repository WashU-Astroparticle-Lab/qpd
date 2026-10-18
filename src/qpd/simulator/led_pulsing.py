"""DAQ-synchronised LED pulses, causal tunneling bursts, and charge resets.

Every call samples one independent chunk. Pulse times are measured from its
start, so resetting the chunk does not change the LED/DAQ phase. All times are
in seconds; ``expected_n_qp`` counts tunneling events, not physical particles.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .quasiparticle_bursts import (
    BurstTruth,
    QuasiparticleBurstModel,
    poisson_burst_times,
)


@dataclass
class LEDTruth:
    """Realised LED pulses and absolute offset-charge resets for one chunk.

    ``bursts`` has one entry per rising edge, including zero-count pulses.
    Pulse end times and tunneling events may extend beyond the chunk: counts
    describe a complete pulse, while the simulator samples only the recording
    window. No pulse or event is carried over into a subsequent chunk.

    ``charge_reset_is_led`` identifies rising-edge resets; the other resets
    belong to an independent homogeneous Poisson process. Values are absolute
    redraws of charge, not additive steps. A reset takes effect at its timestamp.
    """

    pulse_times: np.ndarray
    pulse_end_times: np.ndarray
    bursts: list[BurstTruth]
    initial_n_g: float
    charge_reset_times: np.ndarray
    charge_reset_values: np.ndarray
    charge_reset_is_led: np.ndarray

    @property
    def event_times(self) -> np.ndarray:
        """All LED-induced tunneling times, sorted, including unrecorded tails."""
        events = [burst.event_times for burst in self.bursts if burst.n_qp]
        return np.sort(np.concatenate(events)) if events else np.empty(0)

    def evaluate_n_g(self, t: np.ndarray) -> np.ndarray:
        """Evaluate the piecewise-constant charge, including exact reset times."""
        t = np.asarray(t, dtype=float)
        values = np.concatenate(([self.initial_n_g], self.charge_reset_values))
        return values[np.searchsorted(self.charge_reset_times, t, side="right")]


@dataclass
class LEDPulseModel:
    """Generate an LED pulse train and its tunneling and offset-charge truth.

    The default train starts at 50 ms, stays ON for 50 us, then OFF for 300 ms
    (a 300.05 ms period). Each pulse has a Poisson-distributed tunneling count
    with mean ``expected_n_qp``. Each event's excitation time is uniform across
    the ON window, followed by the existing exponentially modified Gaussian
    delay with defaults mu=1.2 ms, sigma=0.4 ms, and tau=3.7 ms. Negative EMG
    delays are resampled, conditioning the delay distribution to be causal.

    At every rising edge the offset charge is redrawn uniformly over
    ``[n_g_min, n_g_max)``. Independent random resets use the same redraw rule
    at ``charge_jump_rate_hz`` (default one per minute), and charge remains
    fixed between resets. Initial charge is drawn from that same distribution.

    A pulse is included when its rising edge is in ``[0, duration)``. Its full
    ON window and all tunneling events are retained even if they extend past
    the chunk boundary. This model does not include the quiescent parity-flip
    contribution; configure that independently in ``VNASimulator``.
    """

    on_time: float = 50e-6
    off_time: float = 0.300
    first_pulse_time: float = 0.050
    expected_n_qp: float = 10.0
    tau: float = 3.7e-3
    mu: float = 1.2e-3
    sigma: float = 0.4e-3
    charge_jump_rate_hz: float = 1.0 / 60.0
    n_g_min: float = -0.5
    n_g_max: float = 0.5

    def __post_init__(self) -> None:
        for name in (
            "on_time", "off_time", "first_pulse_time", "expected_n_qp",
            "tau", "mu", "sigma", "charge_jump_rate_hz", "n_g_min", "n_g_max",
        ):
            value = np.asarray(getattr(self, name), dtype=float)
            if value.ndim != 0 or not np.isfinite(value):
                raise ValueError(f"{name} must be a finite scalar")
            setattr(self, name, float(value))
        if self.on_time <= 0:
            raise ValueError("on_time must be positive")
        for name in (
            "off_time", "first_pulse_time", "expected_n_qp", "tau", "mu",
            "sigma", "charge_jump_rate_hz",
        ):
            if getattr(self, name) < 0:
                raise ValueError(f"{name} must be non-negative")
        if not np.isfinite(self.period):
            raise ValueError("on_time + off_time must be finite")
        if self.n_g_max <= self.n_g_min or not np.isfinite(
            self.n_g_max - self.n_g_min
        ):
            raise ValueError("n_g_max - n_g_min must be positive and finite")

    @property
    def period(self) -> float:
        """Time between rising edges: ON duration plus OFF duration, seconds."""
        return self.on_time + self.off_time

    def sample(self, duration: float, rng: np.random.Generator) -> LEDTruth:
        """Sample one independent chunk with fixed LED/DAQ alignment.

        Separate random streams keep charge resets unchanged when the pulse
        mean or tunneling profile changes with a fixed upstream RNG seed.
        Zero duration is valid and produces no pulses or reset events.
        """
        duration_value = np.asarray(duration, dtype=float)
        if (
            duration_value.ndim != 0
            or not np.isfinite(duration_value)
            or duration_value < 0
        ):
            raise ValueError("duration must be a finite, non-negative scalar")
        duration = float(duration_value)
        seeds = rng.integers(0, np.iinfo(np.uint64).max, size=4, dtype=np.uint64)
        burst_rng, excitation_rng, reset_time_rng, reset_value_rng = (
            np.random.default_rng(seed) for seed in seeds
        )

        if self.first_pulse_time >= duration:
            pulse_times = np.empty(0, dtype=float)
        else:
            count = int(np.ceil((duration - self.first_pulse_time) / self.period))
            pulse_times = self.first_pulse_time + self.period * np.arange(count)
            pulse_times = pulse_times[pulse_times < duration]

        # Sample at zero before assigning onset times so rejection operates on
        # true delays, without subtracting large absolute timestamps.
        _, raw_bursts = QuasiparticleBurstModel(
            times=np.zeros(pulse_times.size),
            tau=self.tau,
            mu=self.mu,
            sigma=self.sigma,
            expected_n_qp=self.expected_n_qp,
        ).sample(burst_rng)
        bursts: list[BurstTruth] = []
        for onset, burst in zip(pulse_times, raw_bursts):
            delays = burst.event_times.copy()
            negative = delays < 0
            # Nonnegative mu guarantees >= 50% acceptance, including the
            # tau=0 half-normal limit. The existing burst model is unchanged.
            while np.any(negative):
                count = int(np.count_nonzero(negative))
                delays[negative] = burst_rng.normal(
                    self.mu, self.sigma, count
                ) + burst_rng.exponential(self.tau, count)
                negative = delays < 0
            events = np.sort(
                onset + excitation_rng.uniform(0.0, self.on_time, burst.n_qp)
                + delays
            )
            bursts.append(BurstTruth(
                t_arrival=float(onset),
                n_qp=burst.n_qp,
                t_start=float(events[0]) if events.size else float(onset),
                t_end=float(events[-1]) if events.size else float(onset),
                event_times=events,
            ))

        random_reset_times = poisson_burst_times(
            self.charge_jump_rate_hz, duration, reset_time_rng
        )
        reset_times = np.concatenate((pulse_times, random_reset_times))
        reset_is_led = np.concatenate((
            np.ones(pulse_times.size, dtype=bool),
            np.zeros(random_reset_times.size, dtype=bool),
        ))
        order = np.argsort(reset_times, kind="stable")
        initial_n_g = float(reset_value_rng.uniform(self.n_g_min, self.n_g_max))
        reset_values = reset_value_rng.uniform(
            self.n_g_min, self.n_g_max, reset_times.size
        )
        return LEDTruth(
            pulse_times=pulse_times,
            pulse_end_times=pulse_times + self.on_time,
            bursts=bursts,
            initial_n_g=initial_n_g,
            charge_reset_times=reset_times[order],
            charge_reset_values=reset_values[order],
            charge_reset_is_led=reset_is_led[order],
        )


def led_count_sweep(points_per_decade: int = 10) -> np.ndarray:
    """Poisson means: zero plus a log-spaced sweep from 1 through 1000.

    The default gives 31 positive means (ten intervals per decade), plus zero.
    Means are floating point; only the realised per-pulse counts are integers.
    """
    if (
        isinstance(points_per_decade, (bool, np.bool_))
        or not isinstance(points_per_decade, (int, np.integer))
        or points_per_decade <= 0
    ):
        raise ValueError("points_per_decade must be a positive integer")
    return np.concatenate(([0.0], np.geomspace(1.0, 1000.0, 3 * points_per_decade + 1)))
