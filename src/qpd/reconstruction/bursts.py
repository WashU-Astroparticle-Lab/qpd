"""Finding quasiparticle bursts in a reconstructed flip train.

The reconstruction returns tunnelling times, not events. A burst -- an energy
deposition that injects a shower of quasiparticles -- shows up as a *cluster*
of those times, and separating clusters from the Poisson background is a
second, separate inference problem. It is the one that matters for a detector:
the physics question is usually "how many bursts, how big", not "when did flip
1748 happen".

The separation is possible because the two processes have very different
densities. On the reference device the background tunnels at ~10 Hz, so a 12 ms
window holds 0.12 background flips; a 25-quasiparticle burst puts 25 flips in
that same window. Three or four flips inside a few milliseconds is therefore
overwhelming evidence, and the whole detector is a scan statistic saying so
quantitatively.

Two things about this are worth stating up front, because both bite.

* **The detector is blind to what it cannot resolve.** The reconstruction
  cannot separate two flips closer than about a sample, and -- worse -- two
  tunnels in quick succession return the parity to where it started, so they
  are not merely mistimed, they are *invisible*. A burst's recovered
  multiplicity therefore saturates as its true multiplicity grows. That is a
  property of parity readout, not of this clustering; see
  :func:`qpd.reconstruction.benchmark.sweep_burst_size`, which measures it.
* **The background rate is an input, and a wrong one moves the threshold.**
  Take it from the data (``BenchmarkReport.corrected_rate_hz``), not from an
  assumption. Too low a rate manufactures bursts out of background pairs; too
  high a rate dissolves the small ones.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.stats import poisson

__all__ = ["DetectedBurst", "detect_bursts", "match_bursts", "BurstMatch"]


@dataclass
class DetectedBurst:
    """One cluster of flip times judged inconsistent with the background."""

    t_start: float  # [s] first flip in the cluster
    t_end: float  # [s] last flip in the cluster
    n_flips: int  # flips in the cluster -- the recovered multiplicity
    p_value: float  # trials-corrected background probability
    flip_times: np.ndarray = field(
        default_factory=lambda: np.empty(0, float))

    @property
    def span(self) -> float:
        return self.t_end - self.t_start

    @property
    def t_mid(self) -> float:
        return 0.5 * (self.t_start + self.t_end)


def detect_bursts(
    flip_times: np.ndarray,
    background_rate_hz: float,
    *,
    duration: float | None = None,
    max_gap: float | None = None,
    min_flips: int = 3,
    max_p_value: float = 1e-3,
) -> list[DetectedBurst]:
    """Cluster a flip train into bursts, against a Poisson background.

    Flips separated by less than ``max_gap`` are linked into a cluster
    (single-linkage in one dimension), and each cluster is then tested against
    the background: for ``m`` flips spanning ``T`` seconds the probability that
    a Poisson process of rate ``lambda`` puts ``m`` or more events in some
    window of that length is

    ``p = min(1, (D / T) * P[Poisson(lambda * T) >= m])``

    with ``D`` the trace duration. The ``D / T`` factor is the trials
    correction, and it is not optional: the cluster was *selected* for being
    dense, so the uncorrected Poisson tail probability is not a false-alarm
    rate. Without it a 5 s trace of pure background at 10 Hz yields spurious
    "bursts" at a rate of order one per trace.

    Parameters
    ----------
    flip_times : array
        Reconstructed (or true) tunnelling times [s]. Sorted internally.
    background_rate_hz : float
        Poisson rate of the background telegraph. Measure it -- on real data
        use :attr:`~qpd.reconstruction.BenchmarkReport.corrected_rate_hz`,
        which is the decoded rate de-biased by the reconstruction's own
        efficiency and purity.
    duration : float, optional
        Trace duration [s], for the trials factor. Defaults to the span of
        ``flip_times``, which is right for a full trace and slightly
        conservative for a slice.
    max_gap : float, optional
        Linking distance [s]. Defaults to ``0.1 / background_rate_hz`` -- the
        gap the background exceeds 90% of the time, so accidental links are
        rare while genuine intra-burst gaps (sub-millisecond on the reference
        device) are far below it.
    min_flips : int
        Smallest cluster that may be called a burst. Two flips are never
        enough: a background pair straddling ``max_gap`` is common, and no
        p-value on two events survives the trials factor anyway.
    max_p_value : float
        Significance gate after the trials correction.

    Returns
    -------
    list of DetectedBurst, ordered in time.
    """
    t = np.sort(np.asarray(flip_times, dtype=float))
    t = t[np.isfinite(t)]
    rate = float(background_rate_hz)
    if rate <= 0 or not np.isfinite(rate):
        raise ValueError(
            f"background_rate_hz must be positive and finite; got {rate!r}. "
            "Measure it from the trace rather than passing zero.")
    if int(min_flips) < 2:
        raise ValueError("min_flips must be at least 2")
    if t.size < int(min_flips):
        return []

    gap = float(0.1 / rate if max_gap is None else max_gap)
    if gap <= 0:
        raise ValueError("max_gap must be positive")
    span_all = float(t[-1] - t[0])
    dur = float(span_all if duration is None else duration)
    if dur <= 0:
        return []

    # Single-linkage clustering: cut wherever the gap exceeds `gap`.
    cuts = np.flatnonzero(np.diff(t) > gap) + 1
    clusters = np.split(t, cuts) if cuts.size else [t]

    out: list[DetectedBurst] = []
    for c in clusters:
        m = int(c.size)
        if m < int(min_flips):
            continue
        span = float(c[-1] - c[0])
        # A cluster of coincident times has no window to speak of; give it the
        # linking distance so the trials factor stays finite and conservative.
        t_eff = max(span, gap)
        lam = rate * t_eff
        # sf(m - 1) = P[X >= m]
        p_local = float(poisson.sf(m - 1, lam))
        trials = max(dur / t_eff, 1.0)
        p = float(min(1.0, trials * p_local))
        if p <= float(max_p_value):
            out.append(DetectedBurst(
                t_start=float(c[0]), t_end=float(c[-1]), n_flips=m,
                p_value=p, flip_times=c))
    return out


@dataclass
class BurstMatch:
    """One truth burst and the detection (if any) that found it."""

    n_qp_true: int  # quasiparticles the burst actually injected
    n_qp_detected: int  # flips in the matched cluster; 0 if undetected
    t_start: float
    t_end: float
    detected: bool
    p_value: float = float("nan")

    @property
    def bias(self) -> int:
        """Detected minus true multiplicity (negative = flips lost)."""
        return self.n_qp_detected - self.n_qp_true


def match_bursts(
    truth,
    detected: list[DetectedBurst],
    *,
    pad: float = 2e-3,
) -> list[BurstMatch]:
    """Match detected clusters to truth bursts, one to one.

    ``truth`` is the list of :class:`~qpd.simulator.BurstTruth` records the
    simulation produced. A detection matches a truth burst when their intervals
    overlap after padding by ``pad``; where several detections overlap the same
    burst the largest is taken, and each detection is used at most once, so a
    burst split into two clusters counts as found once rather than twice.

    Truth bursts that injected no quasiparticles at all (``n_qp = 0``, which
    Poisson sampling produces for small means) are dropped: there is nothing
    there to detect, and counting them would report an efficiency for events
    that never happened.

    Returns one :class:`BurstMatch` per surviving truth burst.
    """
    dets = sorted(detected, key=lambda d: -d.n_flips)
    used = set()
    out: list[BurstMatch] = []
    for b in truth:
        if int(b.n_qp) <= 0:
            continue
        lo, hi = float(b.t_start) - pad, float(b.t_end) + pad
        hit = None
        for k, d in enumerate(dets):
            if k in used:
                continue
            if d.t_end >= lo and d.t_start <= hi:
                hit = (k, d)
                break
        if hit is None:
            out.append(BurstMatch(
                n_qp_true=int(b.n_qp), n_qp_detected=0,
                t_start=float(b.t_start), t_end=float(b.t_end),
                detected=False))
        else:
            used.add(hit[0])
            out.append(BurstMatch(
                n_qp_true=int(b.n_qp), n_qp_detected=int(hit[1].n_flips),
                t_start=float(b.t_start), t_end=float(b.t_end),
                detected=True, p_value=hit[1].p_value))
    return out
