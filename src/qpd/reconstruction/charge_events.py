"""Offset-charge event detection at a fixed bias point.

At constant ``n_g`` the two parity branches are two stationary blobs, and an
offset-charge jump is the one thing that moves them: both centres slide along
the S21 resonance curve to a new pair of positions, and -- at the reference
geometry, where the drive sits well off the dressed line -- the dominant
observable change is the *separation* (the contrast), with the pair midpoint
nearly fixed. The blobs stay on the projection axis: the perpendicular sagitta
of the move is under 0.15 sigma for any realistic jump, which is why detection
here is one-dimensional and there is no perpendicular-step statistic.

Detection is two-stage: a cheap **screen** proposes candidate boundaries, and
a local HMM likelihood ratio **verifies** each one. Both stages exist because
of measured failures of the single-stage alternatives:

* **The screen cannot compare against the global fit alone.** A global fit
  across a jump does not stay anchored to the pre-jump blobs -- it smears
  toward a compromise, and for a mid-trace jump the misfit against it is
  elevated near-*symmetrically*: a step detector sees a plateau with no step.
  So the screen is a split scan: at each candidate boundary, how much
  nearest-centre distortion ``sum min((x-a)^2, (x-b)^2)`` do freshly fitted
  side pairs remove, compared with the best single pair for the interval? It
  peaks at the jump no matter how the single fit distributes its error.
* **The screen alone cannot be trusted, because refitted statistics are
  occupancy-confounded.** Branch occupancy is correlated over whole dwells,
  so on segment scales it fluctuates wildly (std ~0.17 over half a reference
  trace), and at moderate contrast *any* per-sample statistic whose fits
  adapt per side inherits that: a likelihood gain with the weight pinned at
  1/2 still fakes hundreds of nats on a clean telegraph (the free means slide
  toward the crowded branch), and even nearest-centre distortion fakes
  ~100-500 because the fitted-pair distortion rate itself depends on the
  branch mix. No threshold separates those from a real jump.
* **The verifier prices occupancy correctly by construction.** For each
  candidate, a window of ``min_segment_samples`` on each side is scored by a
  two-state HMM twice -- emission means switching at the boundary (each side's
  own pair) versus one pair for the whole window -- and the gain is the
  difference of *marginal* likelihoods. Under the telegraph prior a lopsided
  stretch of dwells is expected, not evidence, so occupancy-faked candidates
  score ~0; only a genuine change of the emission geometry survives. The
  ratio is insensitive to the exact ``p_flip`` used, since both hypotheses
  share the transition prior.

The distortion and the emissions are normalized by the *minor-axis* sigma of
the global fit, not its EM width: the EM sigma of a jumped trace is inflated
by the very smear being tested for, deflating the statistic exactly when it
matters, while the minor-axis width is orthogonal to blob motion and stays
clean. The acceptance threshold ``min_gain_nats`` applies to the verifier's
likelihood ratio (2 extra fitted means, plus a trials term for the screen)
and is validated against the measured jump-free null in
``checks/study_charge_event_static.py`` rather than derived from an
independence assumption the telegraph does not satisfy.

One case is irreducible: ``delta = +/-0.5`` (mod 1). ``chi_odd(n_g) =
chi_even(n_g + 1/2)`` exactly, so the blob *pair* maps onto itself and the
trace distribution is unchanged; the only effect is one effective branch
relabel at the jump time, indistinguishable from a parity flip in any single
trace. This is the same limit the ramped path documents for its
``delta = +/-0.5`` case, and it costs at most one flip error.

Boundaries may land anywhere down to ~32 samples from the trace ends: the
verifier handles an asymmetric window (a short side still gets its own
HMM-fitted pair), so edge jumps need no separate detector. What a short
segment cannot support is a trustworthy *refit for decoding* -- that rule
(shorter than ``min_segment_samples`` -> declared dead time, never a global
fallback) belongs to the caller,
:func:`~qpd.reconstruction.static_bias.reconstruct_parity_flips_static`,
and is applied uniformly to short end segments and short between-jump
segments alike.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .hmm import forward_backward, gaussian_log_emissions
from .static_bias import StaticBlobModel

__all__ = ["ChargeEventDiagnostics", "detect_charge_events"]

# The screen runs on a decimated trace capped near this many samples, so its
# cost is bounded regardless of trace length; localization and verification
# always use the full-resolution data.
_SCREEN_SAMPLES = 65536
# Cap the candidate grid per interval; localization below the resulting
# stride is recovered by the misfit refinement step.
_MAX_CANDIDATES = 64
# Screening floor (distortion units) and how many of an interval's local
# maxima get a verification attempt. The floor only has to reject obvious
# noise: the verifier is the real threshold. The try budget is generous
# because the screen is occupancy-confounded by design -- fake maxima of
# 100-200 distortion units routinely outrank a true edge jump, and each
# verification costs only ~10 ms, so trying every well-spaced maximum is
# cheaper than trusting the ranking.
_SCREEN_FLOOR = 25.0
_SCREEN_TRIES = 16
# Closest a candidate boundary may sit to an interval end. A foreign stretch
# shorter than this is undetectable against the noise and harmless to the
# decode.
_MIN_SPLIT_MARGIN = 32


@dataclass
class ChargeEventDiagnostics:
    """Everything :func:`detect_charge_events` decided and measured."""

    boundaries: np.ndarray  # sample indices where a new segment starts
    gain_nats: np.ndarray  # verified HMM likelihood-ratio per boundary
    localization_sigma: np.ndarray  # ~1-sigma boundary uncertainty [samples]
    scan: dict = field(default_factory=dict)


def _fit_pair(x: np.ndarray, sigma: float,
              mu_a: float | None = None, mu_b: float | None = None,
              n_iter: int = 25, tol: float = 1e-10
              ) -> tuple[float, float, float]:
    """Two centres minimizing nearest-centre distortion; returns ``-J/2``.

    Lloyd's iteration (2-means in 1-D) from a median-split start unless warm
    centres are given. The score is minus half the distortion in sigma^2
    units so screening gains read like log-likelihood gains.
    """
    if mu_a is None or mu_b is None:
        med = float(np.median(x))
        hi = x >= med
        if not hi.any() or hi.all():
            mu_a = mu_b = float(np.mean(x))
        else:
            mu_a = float(np.mean(x[hi]))
            mu_b = float(np.mean(x[~hi]))
    for _ in range(n_iter):
        use_a = np.abs(x - mu_a) <= np.abs(x - mu_b)
        na = int(np.count_nonzero(use_a))
        if na == 0 or na == x.size:
            # One centre lost its samples: a single blob. Collapse both onto
            # the mean instead of chasing noise splinters.
            mu_a = mu_b = float(np.mean(x))
            break
        new_a = float(np.mean(x[use_a]))
        new_b = float(np.mean(x[~use_a]))
        done = abs(new_a - mu_a) < tol and abs(new_b - mu_b) < tol
        mu_a, mu_b = new_a, new_b
        if done:
            break
    j = float(np.sum(np.minimum((x - mu_a) ** 2, (x - mu_b) ** 2)))
    return mu_a, mu_b, -0.5 * j / (sigma * sigma)


def _pair_nll(x: np.ndarray, mu_a: float, mu_b: float,
              sigma: float) -> np.ndarray:
    """Per-sample negative log-likelihood under a fixed equal-weight pair.

    Parity-blind (both means are allowed) so parity flips and burst-rate
    toggling leave it flat, and occupancy-independent because the pair is
    fixed, not refitted. Crucially this is a *proper scoring rule*, which
    plain nearest-centre distance is not: two spread-out centres score
    better distortion on overlapping data than the true narrow pair (the
    spread catches the tails), so a distance-based comparison between pairs
    of different separations picks the wrong one -- measured as boundary
    refinement sliding to the window edge on a separation-shrink jump. The
    log-density has no such failure. Constant terms (same sigma for both
    pairs) are dropped.
    """
    inv = 0.5 / (sigma * sigma)
    da = (x - mu_a) ** 2
    db = (x - mu_b) ** 2
    m = np.minimum(da, db)
    return m * inv - np.log(0.5 * (1.0 + np.exp(-np.abs(da - db) * inv)))


def _hmm_loglik(x: np.ndarray, mu_a, mu_b, sigma: float,
                p_flip: float) -> float:
    """Marginal log-likelihood of the symmetric two-state chain."""
    _, ll = forward_backward(
        gaussian_log_emissions(x, np.broadcast_to(mu_a, x.shape),
                               np.broadcast_to(mu_b, x.shape), sigma), p_flip)
    return float(ll)


def _hmm_pair(x: np.ndarray, sigma: float, p_flip: float,
              mu_a: float, mu_b: float, n_iter: int = 3
              ) -> tuple[float, float, float]:
    """Branch means fitted by soft-EM under the telegraph prior.

    This -- not a marginal-distribution fit -- is what makes the verifier
    occupancy-proof: the forward-backward posterior assigns whole dwells to
    branches, so a window that happens to sit 80/20 on one branch still gets
    both conditional means right, where a 2-means fit would split the crowded
    blob and orphan the other. Returns the means and the marginal
    log-likelihood evaluated *at* them.
    """
    ll = -np.inf
    for _ in range(n_iter):
        post, ll = forward_backward(
            gaussian_log_emissions(x, np.full(x.size, mu_a),
                                   np.full(x.size, mu_b), sigma), p_flip)
        wa = 1.0 - post
        sa, sb = float(wa.sum()), float(post.sum())
        if sa < 1.0 or sb < 1.0:
            break
        mu_a = float((wa * x).sum() / sa)
        mu_b = float((post * x).sum() / sb)
    return mu_a, mu_b, _hmm_loglik(x, mu_a, mu_b, sigma, p_flip)


def _verify_candidate(x: np.ndarray, sigma: float, b: int, lo: int, hi: int,
                      window: int, p_flip: float,
                      init: tuple[float, float]
                      ) -> tuple[float, tuple[float, float],
                                 tuple[float, float]]:
    """HMM likelihood ratio for a boundary at ``b``: switch-at-b vs one pair.

    Both hypotheses share the transition prior and the same fitting
    procedure; only the emission means differ -- each side's own pair against
    one pair for the whole window. Dwell-correlated occupancy is priced
    identically on both sides of the ratio, so the null gain is just the
    cost of two extra fitted means. Returns the gain and the two side pairs
    (for refinement).
    """
    w_lo = max(lo, b - window)
    w_hi = min(hi, b + window)
    if b - w_lo < 16 or w_hi - b < 16:
        return -np.inf, init, init
    sa, sb, ll_single = _hmm_pair(x[w_lo:w_hi], sigma, p_flip, *init)
    la, lb, ll_l = _hmm_pair(x[w_lo:b], sigma, p_flip, sa, sb)
    ra, rb, ll_r = _hmm_pair(x[b:w_hi], sigma, p_flip, sa, sb)
    # The split model's likelihood is evaluated jointly (one chain, means
    # switching at b) so the two hypotheses integrate over the same paths.
    n_l = b - w_lo
    n_r = w_hi - b
    mu_a = np.concatenate([np.full(n_l, la), np.full(n_r, ra)])
    mu_b = np.concatenate([np.full(n_l, lb), np.full(n_r, rb)])
    ll_split = _hmm_loglik(x[w_lo:w_hi], mu_a, mu_b, sigma, p_flip)
    return ll_split - ll_single, (la, lb), (ra, rb)


def _refine_boundary(x: np.ndarray, sigma: float, lo: int, hi: int,
                     left: tuple[float, float], right: tuple[float, float],
                     ) -> tuple[int, float]:
    """Exact change-point inside ``[lo, hi)`` by cumulative pair likelihood.

    Returns the boundary and an approximate 1-sigma localization width: the
    half-width of the region where the switching cost stays within 2 nats
    of its minimum. Small jumps localize poorly (the cost trough is flat)
    and the width says so.
    """
    xs = x[lo:hi]
    ma = _pair_nll(xs, left[0], left[1], sigma)
    mb = _pair_nll(xs, right[0], right[1], sigma)
    total = np.concatenate(([0.0], np.cumsum(ma))) + (
        np.sum(mb) - np.concatenate(([0.0], np.cumsum(mb)))
    )
    k = int(np.argmin(total))
    width = 0.5 * np.count_nonzero(total <= total[k] + 2.0)
    return lo + k, float(max(width, 0.5))


def _scan_interval(x: np.ndarray, sigma: float, lo: int, hi: int,
                   stride: int, decim: int):
    """Screening gains for every candidate split of ``[lo, hi)``.

    Runs on the trace decimated by ``decim``; side fits are warm-started from
    the previous candidate (the centres drift slowly as the boundary moves).
    Returns ``(ks, gains)`` on the full-resolution index grid. Candidates run
    to within ``_MIN_SPLIT_MARGIN`` of the interval ends: a short side cannot
    be *refit for decoding*, but it can be screened and verified, which is
    how edge jumps are found without any separate machinery.
    """
    ks = np.arange(lo + _MIN_SPLIT_MARGIN, hi - _MIN_SPLIT_MARGIN + 1,
                   stride, dtype=int)
    if ks.size == 0:
        return ks, np.empty(0)
    xs = x[lo:hi:decim]
    _, _, s_whole = _fit_pair(xs, sigma)
    gains = np.empty(ks.size)
    la = lb = ra = rb = None
    for i, k in enumerate(ks):
        kd = max(1, min((k - lo) // decim, xs.size - 1))
        la, lb, s_l = _fit_pair(xs[:kd], sigma, la, lb, n_iter=8)
        ra, rb, s_r = _fit_pair(xs[kd:], sigma, ra, rb, n_iter=8)
        gains[i] = s_l + s_r - s_whole
    return ks, gains


def _local_maxima(gains: np.ndarray) -> np.ndarray:
    """Indices of local maxima, best first."""
    if gains.size == 0:
        return np.empty(0, dtype=int)
    if gains.size == 1:
        return np.array([0])
    left_ok = np.r_[True, gains[1:] >= gains[:-1]]
    right_ok = np.r_[gains[:-1] >= gains[1:], True]
    idx = np.flatnonzero(left_ok & right_ok)
    return idx[np.argsort(gains[idx])[::-1]]


def detect_charge_events(
    iq: np.ndarray,
    model: StaticBlobModel,
    dt: float,
    *,
    scan_stride: int = 512,
    min_segment_samples: int = 2500,
    min_gain_nats: float = 15.0,
    max_events: int = 8,
    p_flip: float = 1e-3,
) -> ChargeEventDiagnostics:
    """Find offset-charge jumps in a fixed-bias trace, blind.

    Parameters
    ----------
    iq : np.ndarray
        Complex trace.
    model : StaticBlobModel
        The global two-blob fit; supplies the projection axis and the
        minor-axis noise width (fallback: the EM width).
    dt : float
        Sample spacing [s]; carried into the diagnostics only.
    scan_stride : int
        Candidate-boundary spacing [samples]. On long traces the effective
        stride grows so the grid stays at most ~64 candidates per interval;
        localization finer than the stride comes from the misfit refinement.
    min_segment_samples : int
        Per-side length of the verification window: enough dwells on each
        branch to pin two conditional means (0.25 s at 10 kSa/s and ~17 Hz
        is about four expected dwells). Boundaries themselves may land much
        closer to the trace ends; whether the resulting short segment can be
        refit for decoding is the caller's rule, not the detector's.
    min_gain_nats : float
        Acceptance threshold on the verifier's HMM likelihood ratio.
        Validated against the measured jump-free null in
        ``checks/study_charge_event_static.py``.
    max_events : int
        Hard cap on accepted boundaries. Real charge events are rare
        (< ~0.1 Hz); hitting the cap is itself a sign that something other
        than charge jumps is moving the blobs, and is recorded in ``scan``.
    p_flip : float
        Per-sample flip probability for the verifier's transition prior. The
        likelihood *ratio* is insensitive to it (both hypotheses share the
        prior), so an order of magnitude is enough.

    Returns
    -------
    ChargeEventDiagnostics
        Accepted boundaries (sample indices, sorted), their verified gains
        and localization widths, and the scan record -- including
        ``max_gain_nats``, the largest *verification* gain seen anywhere
        (accepted or not), which is what the null calibration measures.
    """
    x = model.project(np.asarray(iq))
    n = x.size
    sigma = float(model.diagnostics.get("sigma_minor_axis", model.sigma)
                  or model.sigma)
    min_seg = int(min_segment_samples)
    decim = max(1, int(np.ceil(n / _SCREEN_SAMPLES)))

    accepted: list[tuple[int, float, float]] = []
    max_gain_seen = -np.inf
    max_screen_seen = -np.inf
    n_scans = n_verifications = 0
    intervals = [(0, n)]
    init = (float(model.mu_a), float(model.mu_b))
    while intervals and len(accepted) < int(max_events):
        lo, hi = intervals.pop()
        stride = max(int(scan_stride),
                     int(np.ceil((hi - lo) / _MAX_CANDIDATES)))
        ks, gains = _scan_interval(x, sigma, lo, hi, stride, decim)
        n_scans += 1
        if gains.size:
            max_screen_seen = max(max_screen_seen, float(gains.max()))
        tried: list[int] = []
        for idx in _local_maxima(gains):
            if len(tried) >= _SCREEN_TRIES or gains[idx] <= _SCREEN_FLOOR:
                break
            k = int(ks[idx])
            # Suppress same-bump duplicates only: a coarser spacing (say the
            # verification window) lets a strong occupancy fake shadow a
            # genuine edge candidate sitting within it.
            if any(abs(k - t) <= stride for t in tried):
                continue
            tried.append(k)
            gain, left, right = _verify_candidate(x, sigma, k, lo, hi,
                                                  min_seg, p_flip, init)
            n_verifications += 1
            if gain > max_gain_seen:
                max_gain_seen = gain
            if gain <= float(min_gain_nats):
                continue
            # A candidate near a real jump verifies even when misplaced (the
            # window straddles the jump either way), and one refinement pass
            # under side pairs fitted to *mixed* windows inherits their
            # compromise. Alternate: refine the boundary, refit the side
            # pairs around it, repeat -- an EM-like loop that walks onto the
            # true change point. The radius must cover the whole
            # verification window, or a candidate a couple of strides off
            # walks in short steps and runs out of rounds.
            radius = max(stride, min_seg)
            b, width = k, float(stride)
            for _ in range(6):
                b_new, width = _refine_boundary(x, sigma,
                                                max(lo, b - radius),
                                                min(hi, b + radius),
                                                left, right)
                moved = abs(b_new - b)
                b = b_new
                if moved <= 1:
                    break
                l_lo = max(lo, b - min_seg)
                r_hi = min(hi, b + min_seg)
                if b - l_lo >= 16:
                    left = _hmm_pair(x[l_lo:b], sigma, p_flip, *left)[:2]
                if r_hi - b >= 16:
                    right = _hmm_pair(x[b:r_hi], sigma, p_flip, *right)[:2]
            # Re-verify at the refined boundary: the recorded gain must
            # belong to the boundary actually reported.
            gain, _, _ = _verify_candidate(x, sigma, b, lo, hi, min_seg,
                                           p_flip, init)
            n_verifications += 1
            if gain > max_gain_seen:
                max_gain_seen = gain
            if gain > float(min_gain_nats):
                accepted.append((b, float(gain), width))
                intervals.append((lo, b))
                intervals.append((b, hi))
                break

    # Prune against final neighbours. During the search a boundary is
    # verified inside its interval *at the time*, and a misplaced candidate
    # that straddles a real jump passes that test. Once every boundary is
    # placed, each one must justify itself with the window clipped to its
    # neighbours: next to the correctly placed jump, a spurious straddler's
    # window is jump-free and its gain collapses. Weakest first, repeated
    # until stable; surviving gains are the honest neighbour-conditional
    # ones.
    accepted.sort(key=lambda a: a[0])
    n_pruned = 0
    changed = True
    while changed and accepted:
        changed = False
        bs = [a[0] for a in accepted]
        for i in sorted(range(len(accepted)), key=lambda j: accepted[j][1]):
            lo_i = bs[i - 1] if i > 0 else 0
            hi_i = bs[i + 1] if i + 1 < len(bs) else n
            gain, _, _ = _verify_candidate(x, sigma, bs[i], lo_i, hi_i,
                                           min_seg, p_flip, init)
            n_verifications += 1
            if gain <= float(min_gain_nats):
                accepted.pop(i)
                n_pruned += 1
                changed = True
                break
            accepted[i] = (bs[i], float(gain), accepted[i][2])

    boundaries = np.array([a[0] for a in accepted], dtype=int)
    gains_out = np.array([a[1] for a in accepted], dtype=float)
    widths = np.array([a[2] for a in accepted], dtype=float)

    return ChargeEventDiagnostics(
        boundaries=boundaries,
        gain_nats=gains_out,
        localization_sigma=widths,
        scan={
            "stride": max(int(scan_stride), int(np.ceil(n / _MAX_CANDIDATES))),
            "decimation": decim,
            "n_scans": n_scans,
            "n_verifications": n_verifications,
            "max_gain_nats": float(max_gain_seen),
            "max_screen_gain": float(max_screen_seen),
            "n_pruned": n_pruned,
            "min_gain_nats": float(min_gain_nats),
            "sigma": sigma,
            "hit_max_events": len(accepted) >= int(max_events),
            "dt": float(dt),
        },
    )
