"""Continuous-time two-state Markov parity trajectory (Gillespie)."""

from __future__ import annotations

from typing import Literal

import numpy as np

EVEN = np.int8(0)
ODD = np.int8(1)


def parity_from_flip_times(
    t_grid: np.ndarray,
    flip_times: np.ndarray,
    initial_state: int,
) -> np.ndarray:
    """Map a set of parity-flip times onto a time grid.

    Every entry in `flip_times` toggles parity. The state at each grid point is
    ``initial_state XOR (number of flips with flip_time <= t) mod 2``.

    `flip_times` need not be sorted (sorted internally). Flips at or before
    `t_grid[0]` shift the visible state by a constant; flips after `t_grid[-1]`
    have no effect — both handled without special-casing, so flip times outside
    the grid are fine. `side="right"` (a flip exactly on a grid sample takes
    effect at that sample) matches the convention in
    `offset_charge.ChargeJumpEvents.evaluate`.

    Returns an int8 array of the same length as `t_grid`.
    """
    t_grid = np.asarray(t_grid, dtype=float)
    if t_grid.size == 0:
        return np.empty(0, dtype=np.int8)
    ft = np.asarray(flip_times, dtype=float)
    if ft.size == 0:
        return np.full(t_grid.size, np.int8(initial_state), dtype=np.int8)
    ft = np.sort(ft)
    n_flips_le = np.searchsorted(ft, t_grid, side="right")
    return (np.int8(initial_state) ^ (n_flips_le & 1).astype(np.int8)).astype(
        np.int8
    )


def generate_parity_trajectory(
    t_grid: np.ndarray,
    gamma_e_to_o: float,
    gamma_o_to_e: float,
    rng: np.random.Generator,
    initial: Literal["steady_state", "even", "odd"] = "steady_state",
    extra_flip_times: np.ndarray | None = None,
    return_flip_times: bool = False,
):
    """Sample a two-state CTMC on a uniform time grid.

    Switching events are drawn from exponential waiting times with rate
    `gamma_e_to_o` while in the even state and `gamma_o_to_e` while in the
    odd state, then evaluated on `t_grid` via `searchsorted`.

    `extra_flip_times`, if given, are additional parity-flip times (e.g.
    quasiparticle-burst tunneling events from
    :class:`qpd.simulator.quasiparticle_bursts.QuasiparticleBurstModel`)
    superimposed on the background CTMC: the CTMC switches and the extra flips
    are merged into a single sorted set of toggles. This treats the two as
    independent point processes — the CTMC's next dwell rate is NOT re-derived
    from the post-flip state after a burst event. That is an acceptable
    approximation when bursts are sparse and fast relative to the background
    telegraph.

    Returns an int8 array of the same length as `t_grid` containing 0
    (even) or 1 (odd).

    If `return_flip_times` is True, returns ``(parity, flip_times)`` where
    ``flip_times`` is the sorted float array of every parity-toggle time that
    falls inside the grid window — the merged set of background CTMC switches
    and burst (``extra_flip_times``) events. This is the exact, sub-sample
    ground truth for parity-flip timing (the grid-quantised transitions in the
    returned array would lose flips separated by less than one sample).
    """
    t_grid = np.asarray(t_grid, dtype=float)
    if t_grid.size == 0:
        empty = np.empty(0, dtype=np.int8)
        if return_flip_times:
            return empty, np.empty(0, dtype=float)
        return empty

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

    # Burst (extra) flips are an EXTERNAL forcing interleaved with the
    # background dwell loop: at each burst time the parity toggles and the
    # background CTMC resamples its next dwell from the post-burst state. The
    # exponential's memorylessness makes discarding the in-flight dwell and
    # resampling exact, so the background re-equilibrates to its steady state
    # after a burst rather than being globally inverted (which an odd number of
    # toggles overlaid on a fixed trajectory would do).
    if extra_flip_times is None:
        extra = np.empty(0, dtype=float)
    else:
        extra = np.sort(np.asarray(extra_flip_times, dtype=float))

    # Burst flips at or before the grid start just shift the initial parity;
    # those past the grid end never affect the visible trajectory.
    if extra.size:
        n_pre = int(np.searchsorted(extra, t_start, side="right"))
        if n_pre & 1:
            initial_state = ODD if initial_state == EVEN else EVEN
        extra = extra[(extra > t_start) & (extra <= t_end)]

    state = initial_state
    switch_times: list[float] = []
    t = t_start
    bi = 0
    while True:
        rate = gamma_e_to_o if state == EVEN else gamma_o_to_e
        if rate <= 0:
            # Background frozen; only remaining burst flips still toggle parity.
            while bi < extra.size:
                switch_times.append(float(extra[bi]))
                state = ODD if state == EVEN else EVEN
                bi += 1
            break
        dwell = rng.exponential(1.0 / rate)
        t_bg = t + dwell
        if bi < extra.size and extra[bi] <= t_bg:
            # A burst preempts this dwell: toggle at the burst time, then
            # resample the next background dwell from the new state.
            tb = float(extra[bi])
            bi += 1
            state = ODD if state == EVEN else EVEN
            switch_times.append(tb)
            t = tb
            continue
        if t_bg > t_end:
            break
        state = ODD if state == EVEN else EVEN
        switch_times.append(t_bg)
        t = t_bg

    flips = np.sort(np.asarray(switch_times, dtype=float))
    parity = parity_from_flip_times(t_grid, flips, int(initial_state))
    if return_flip_times:
        return parity, flips
    return parity
