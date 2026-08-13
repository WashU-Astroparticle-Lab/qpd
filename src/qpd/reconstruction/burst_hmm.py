"""Burst-aware parity decoding: a parity x regime hidden Markov model.

The two-state decoder in :mod:`.hmm` fits **one** flip probability to the
whole trace. A quasiparticle burst is a local excursion of the tunnelling
rate by ~3 orders of magnitude, and a single global rate cannot represent
it: the background dwell dominates the fit, Viterbi prices every in-burst
toggle as a rare event, and the recovered multiplicity saturates near the
cluster threshold regardless of the true burst size (issue #40).

This module extends the chain with a latent *regime* bit:

======  ========  ========
state   parity    regime
======  ========  ========
0       A         quiet
1       B         quiet
2       A         burst
3       B         burst
======  ========  ========

Emissions depend on parity alone -- the readout cannot see the regime, only
the branch -- so the fitted blob model is reused unchanged. All the physics
sits in the transition matrix:

* parity flips with probability ``p_quiet`` per sample in the quiet regime
  (the background telegraph) and ``p_burst`` in the burst regime;
* the regime enters a burst with probability ``epsilon`` per sample and
  leaves with ``dt / burst_tau``, the burst model's own decay.

Three of the four parameters are not free: ``p_quiet`` is re-estimated from
the trace exactly as the two-state decoder estimates its global rate (and,
usefully, is no longer inflated by the bursts it now explains separately),
``burst_tau`` is the known quasiparticle decay time, and ``p_burst`` is
**pinned** -- letting EM drive it toward 1/2 would turn the burst regime
into a "don't care" state that absorbs noise, so it is a fixed operating
point, not a fit.

``epsilon`` deserves a word, because it looks like an assumed burst rate
and is not. It enters the decode once per burst as an entry cost of
``ln(1/epsilon)`` nats, against evidence of ~``ln(p_burst/p_quiet)`` nats
*per recovered flip* -- so an order-of-magnitude error in ``epsilon`` moves
the detection threshold by a fraction of a flip. It is a false-alarm
threshold in the same sense as the scan statistic's p-value gate in
:mod:`.bursts`, which remains the reported significance either way.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from math import log

import numpy as np

from .hmm import _MIN_P, decode_with_rate, gaussian_log_emissions

__all__ = [
    "BurstAwareResult",
    "burst_transition_matrix",
    "forward_backward_regime",
    "viterbi_regime",
    "decode_burst_aware",
]

# State layout: parity bit and regime bit of each of the four chain states.
_PARITY = np.array([0, 1, 0, 1])
_REGIME = np.array([0, 0, 1, 1])


@dataclass
class BurstAwareResult:
    """Posterior and MAP decoding of the parity x regime chain."""

    parity_posterior: np.ndarray  # (n,) P[parity = B | entire trace]
    burst_posterior: np.ndarray  # (n,) P[regime = burst | entire trace]
    parity_path: np.ndarray  # (n,) int8 parity bit of the Viterbi path
    regime_path: np.ndarray  # (n,) int8 regime bit of the Viterbi path
    log_likelihood: float
    p_quiet: float  # fitted per-sample flip probability, quiet regime
    p_burst: float  # pinned per-sample flip probability, burst regime
    diagnostics: dict = field(default_factory=dict)


def burst_transition_matrix(
    p_quiet: float,
    p_burst: float,
    entry: float,
    exit_: float,
) -> np.ndarray:
    """Per-sample transition matrix of the four-state chain.

    The regime move and the parity move factorise within one step, with the
    flip probability taken from the *source* regime -- a flip coincident
    with regime entry is charged at the quiet rate, which costs at most one
    flip of evidence at the burst edge and keeps the factorisation exact.
    """
    pq = float(min(max(p_quiet, _MIN_P), 0.5))
    # A burst regime slower than the background would invert the model's
    # meaning; clip rather than fail so a pathological fit degrades softly.
    pb = float(min(max(p_burst, pq), 0.5))
    en = float(min(max(entry, _MIN_P), 0.5))
    ex = float(min(max(exit_, _MIN_P), 0.5))
    regime = np.array([[1.0 - en, en], [ex, 1.0 - ex]])
    t = np.empty((4, 4))
    for i in range(4):
        p = pb if _REGIME[i] else pq
        for j in range(4):
            flip = p if _PARITY[i] != _PARITY[j] else 1.0 - p
            t[i, j] = regime[_REGIME[i], _REGIME[j]] * flip
    return t


def _scaled_emissions4(log_emit: np.ndarray) -> tuple[np.ndarray, float]:
    """Four-row emission likelihoods rescaled per sample against underflow."""
    peak = log_emit.max(axis=0)
    return np.exp(log_emit - peak), float(peak.sum())


def forward_backward_regime(
    log_emit4: np.ndarray, trans: np.ndarray
) -> tuple[np.ndarray, float]:
    """Posterior over the four states and the total log-likelihood.

    Same scaled-domain recursion as the two-state
    :func:`~qpd.reconstruction.hmm.forward_backward`, specialised to scalar
    arithmetic for the same reason: per-step NumPy calls cost an order of
    magnitude more than the sixteen multiplies they replace.

    Returns ``(posterior, log_likelihood)`` with posterior of shape
    ``(4, n)``.
    """
    e, log_scale = _scaled_emissions4(log_emit4)
    e0, e1, e2, e3 = (row.tolist() for row in e)
    n = len(e0)
    (t00, t01, t02, t03), (t10, t11, t12, t13), \
        (t20, t21, t22, t23), (t30, t31, t32, t33) = np.asarray(
            trans, dtype=float).tolist()

    fwd0 = [0.0] * n
    fwd1 = [0.0] * n
    fwd2 = [0.0] * n
    fwd3 = [0.0] * n
    a0, a1, a2, a3 = 0.25 * e0[0], 0.25 * e1[0], 0.25 * e2[0], 0.25 * e3[0]
    s = a0 + a1 + a2 + a3
    log_scale += log(s)
    a0 /= s
    a1 /= s
    a2 /= s
    a3 /= s
    fwd0[0], fwd1[0], fwd2[0], fwd3[0] = a0, a1, a2, a3
    for k in range(1, n):
        b0 = (a0 * t00 + a1 * t10 + a2 * t20 + a3 * t30) * e0[k]
        b1 = (a0 * t01 + a1 * t11 + a2 * t21 + a3 * t31) * e1[k]
        b2 = (a0 * t02 + a1 * t12 + a2 * t22 + a3 * t32) * e2[k]
        b3 = (a0 * t03 + a1 * t13 + a2 * t23 + a3 * t33) * e3[k]
        s = b0 + b1 + b2 + b3
        log_scale += log(s)
        a0, a1, a2, a3 = b0 / s, b1 / s, b2 / s, b3 / s
        fwd0[k], fwd1[k], fwd2[k], fwd3[k] = a0, a1, a2, a3

    post = np.empty((4, n))
    b0 = b1 = b2 = b3 = 1.0
    post[0, n - 1] = fwd0[n - 1]
    post[1, n - 1] = fwd1[n - 1]
    post[2, n - 1] = fwd2[n - 1]
    post[3, n - 1] = fwd3[n - 1]
    for k in range(n - 2, -1, -1):
        c0 = b0 * e0[k + 1]
        c1 = b1 * e1[k + 1]
        c2 = b2 * e2[k + 1]
        c3 = b3 * e3[k + 1]
        b0 = t00 * c0 + t01 * c1 + t02 * c2 + t03 * c3
        b1 = t10 * c0 + t11 * c1 + t12 * c2 + t13 * c3
        b2 = t20 * c0 + t21 * c1 + t22 * c2 + t23 * c3
        b3 = t30 * c0 + t31 * c1 + t32 * c2 + t33 * c3
        s = b0 + b1 + b2 + b3
        b0 /= s
        b1 /= s
        b2 /= s
        b3 /= s
        g0 = fwd0[k] * b0
        g1 = fwd1[k] * b1
        g2 = fwd2[k] * b2
        g3 = fwd3[k] * b3
        s = g0 + g1 + g2 + g3
        post[0, k] = g0 / s
        post[1, k] = g1 / s
        post[2, k] = g2 / s
        post[3, k] = g3 / s
    return post, log_scale


def viterbi_regime(log_emit4: np.ndarray, trans: np.ndarray) -> np.ndarray:
    """Most likely four-state sequence, as an int8 array of state indices.

    Scalar-specialised like the forward-backward above; the sixteen
    max-comparisons per step beat a NumPy argmax on a (4, 4) array by an
    order of magnitude.
    """
    lt = np.log(np.asarray(trans, dtype=float)).tolist()
    (t00, t01, t02, t03), (t10, t11, t12, t13), \
        (t20, t21, t22, t23), (t30, t31, t32, t33) = lt
    l0, l1, l2, l3 = (row.tolist() for row in log_emit4)
    n = len(l0)

    back = [[0] * 4 for _ in range(n)]
    d0 = log(0.25) + l0[0]
    d1 = log(0.25) + l1[0]
    d2 = log(0.25) + l2[0]
    d3 = log(0.25) + l3[0]
    for k in range(1, n):
        bk = back[k]

        m, i = d0 + t00, 0
        v = d1 + t10
        if v > m:
            m, i = v, 1
        v = d2 + t20
        if v > m:
            m, i = v, 2
        v = d3 + t30
        if v > m:
            m, i = v, 3
        n0, bk[0] = m + l0[k], i

        m, i = d0 + t01, 0
        v = d1 + t11
        if v > m:
            m, i = v, 1
        v = d2 + t21
        if v > m:
            m, i = v, 2
        v = d3 + t31
        if v > m:
            m, i = v, 3
        n1, bk[1] = m + l1[k], i

        m, i = d0 + t02, 0
        v = d1 + t12
        if v > m:
            m, i = v, 1
        v = d2 + t22
        if v > m:
            m, i = v, 2
        v = d3 + t32
        if v > m:
            m, i = v, 3
        n2, bk[2] = m + l2[k], i

        m, i = d0 + t03, 0
        v = d1 + t13
        if v > m:
            m, i = v, 1
        v = d2 + t23
        if v > m:
            m, i = v, 2
        v = d3 + t33
        if v > m:
            m, i = v, 3
        n3, bk[3] = m + l3[k], i

        d0, d1, d2, d3 = n0, n1, n2, n3

    path = np.empty(n, dtype=np.int8)
    state = max(range(4), key=(d0, d1, d2, d3).__getitem__)
    path[n - 1] = state
    for k in range(n - 1, 0, -1):
        state = back[k][state]
        path[k - 1] = state
    return path


def decode_burst_aware(
    x: np.ndarray,
    mu_a: np.ndarray,
    mu_b: np.ndarray,
    sigma: float,
    dt: float,
    p_flip_init: float = 1e-3,
    n_iter: int = 4,
    p_burst: float = 0.3,
    burst_rate_hz: float = 1.0,
    burst_tau: float = 1e-3,
) -> BurstAwareResult:
    """Decode a static-bias trace under the parity x regime model.

    The quiet-regime flip probability is estimated from the data by the same
    hard-EM loop the two-state decoder uses, except the counting is restricted
    to samples the current Viterbi path assigns to the quiet regime -- so
    bursts stop inflating the background rate estimate as they do in the
    global fit. The two-state estimate seeds the iteration.

    Parameters
    ----------
    x, mu_a, mu_b, sigma
        Projected trace and per-sample branch means, as in
        :func:`~qpd.reconstruction.hmm.decode_with_rate`.
    dt : float
        Sample period [s]; the regime rates below are physical and need it.
    p_flip_init : float
        Seed for the quiet-rate estimate.
    n_iter : int
        Hard-EM refinements of ``p_quiet``.
    p_burst : float
        Pinned per-sample flip probability inside a burst. Not fitted -- see
        the module docstring for why.
    burst_rate_hz : float
        Expected burst occurrence rate [Hz]; sets the regime entry probability
        ``epsilon = burst_rate_hz * dt``. Logarithmic sensitivity: a factor of
        100 here moves the detection threshold by about one flip of evidence.
    burst_tau : float
        Quasiparticle burst decay time [s]; sets the regime exit probability
        ``dt / burst_tau``.
    """
    x = np.asarray(x)
    log_emit2 = gaussian_log_emissions(x, mu_a, mu_b, sigma)
    log_emit4 = log_emit2[[0, 1, 0, 1]]
    entry = float(burst_rate_hz) * dt
    exit_ = dt / float(burst_tau)

    # Seed the quiet rate with the global two-state estimate. It is biased
    # high when bursts are present, which is the safe direction to start from:
    # the first regime decode then reassigns the bursts and the estimate
    # relaxes to the true background.
    _, p_quiet, seed_history = decode_with_rate(
        x, mu_a, mu_b, sigma, p_flip_init=p_flip_init, n_iter=max(1, n_iter))

    history = [p_quiet]
    path4 = None
    for _ in range(max(1, int(n_iter))):
        trans = burst_transition_matrix(p_quiet, p_burst, entry, exit_)
        path4 = viterbi_regime(log_emit4, trans)
        quiet = _REGIME[path4[:-1]] == 0
        n_quiet = int(np.count_nonzero(quiet))
        if n_quiet == 0:
            break
        flips = np.count_nonzero(
            (_PARITY[path4[1:]] != _PARITY[path4[:-1]]) & quiet)
        p_new = max(flips / n_quiet, 1e-9)
        history.append(p_new)
        converged = abs(p_new - p_quiet) < 1e-3 * p_quiet
        p_quiet = p_new
        if converged:
            break

    trans = burst_transition_matrix(p_quiet, p_burst, entry, exit_)
    post4, loglik = forward_backward_regime(log_emit4, trans)
    path4 = viterbi_regime(log_emit4, trans)

    return BurstAwareResult(
        parity_posterior=post4[1] + post4[3],
        burst_posterior=post4[2] + post4[3],
        parity_path=_PARITY[path4].astype(np.int8),
        regime_path=_REGIME[path4].astype(np.int8),
        log_likelihood=loglik,
        p_quiet=float(p_quiet),
        p_burst=float(p_burst),
        diagnostics={
            "epsilon": entry,
            "exit_prob": exit_,
            "burst_rate_hz": float(burst_rate_hz),
            "burst_tau": float(burst_tau),
            "p_quiet_history": history,
            "p_global_seed": seed_history[-1],
            "state_path": path4,
        },
    )
