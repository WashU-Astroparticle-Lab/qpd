"""Continuous-time two-state Markov parity trajectory (Gillespie)."""

from __future__ import annotations

from typing import Literal

import numpy as np

EVEN = np.int8(0)
ODD = np.int8(1)


def generate_parity_trajectory(
    t_grid: np.ndarray,
    gamma_e_to_o: float,
    gamma_o_to_e: float,
    rng: np.random.Generator,
    initial: Literal["steady_state", "even", "odd"] = "steady_state",
) -> np.ndarray:
    """Sample a two-state CTMC on a uniform time grid.

    Switching events are drawn from exponential waiting times with rate
    `gamma_e_to_o` while in the even state and `gamma_o_to_e` while in the
    odd state, then evaluated on `t_grid` via `searchsorted`.

    Returns an int8 array of the same length as `t_grid` containing 0
    (even) or 1 (odd).
    """
    t_grid = np.asarray(t_grid, dtype=float)
    if t_grid.size == 0:
        return np.empty(0, dtype=np.int8)

    if gamma_e_to_o < 0 or gamma_o_to_e < 0:
        raise ValueError("rates must be non-negative")

    t_start = float(t_grid[0])
    t_end = float(t_grid[-1])

    if initial == "even":
        initial_state = EVEN
    elif initial == "odd":
        initial_state = ODD
    elif initial == "steady_state":
        total = gamma_e_to_o + gamma_o_to_e
        p_even = gamma_o_to_e / total if total > 0 else 0.5
        initial_state = EVEN if rng.random() < p_even else ODD
    else:
        raise ValueError(f"unknown initial: {initial!r}")

    state = initial_state
    switch_times: list[float] = []
    switch_states: list[int] = []
    t = t_start
    while True:
        rate = gamma_e_to_o if state == EVEN else gamma_o_to_e
        if rate <= 0:
            break
        dwell = rng.exponential(1.0 / rate)
        t = t + dwell
        if t > t_end:
            break
        state = ODD if state == EVEN else EVEN
        switch_times.append(t)
        switch_states.append(int(state))

    if not switch_times:
        return np.full(t_grid.size, initial_state, dtype=np.int8)

    # Interval i covers [boundaries[i], boundaries[i+1]); state in that
    # interval is states_seq[i]. boundary[0] is a sentinel before t_start.
    states_seq = np.concatenate(
        ([initial_state], np.asarray(switch_states, dtype=np.int8))
    )
    boundaries = np.concatenate(([t_start - 1.0], np.asarray(switch_times)))
    idx = np.searchsorted(boundaries, t_grid, side="right") - 1
    return states_seq[idx]
