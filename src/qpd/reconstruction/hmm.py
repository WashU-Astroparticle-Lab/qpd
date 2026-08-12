"""Two-state HMM inference with time-varying Gaussian emissions.

The parity trajectory is a two-state telegraph process observed through
additive complex Gaussian noise. Both ingredients of the HMM are supplied
by the caller:

* the **transition model** -- per-sample flip probability ``p = Gamma*dt``,
  which encodes "flips are rare on the sample timescale" and is what keeps
  single noisy samples from being read as flips;
* the **emission model** -- for every sample, the two candidate means
  ``mu_A(t)``, ``mu_B(t)`` that the signal would sit near under each branch
  (these move with n_g), plus the noise sigma.

:func:`forward_backward` returns the per-sample posterior over branches
given the *whole* trace, and :func:`viterbi` returns the single most likely
branch sequence. Together these integrate evidence across samples weighted
by the dwell-time prior -- effectively a matched filter whose averaging
length adapts to the local branch separation, which is what makes
per-sample SNR of order unity workable.

Both recursions are exact and O(n). The two-state case is specialised to
scalar arithmetic (a length-n NumPy loop would be an order of magnitude
slower); :func:`forward_backward_reference` is a transparent general-K
log-domain implementation kept for testing the fast path.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from math import log

import numpy as np

__all__ = [
    "HMMResult",
    "gaussian_log_emissions",
    "forward_backward",
    "viterbi",
    "decode",
    "decode_with_rate",
    "decoded_path",
    "forward_backward_reference",
]

DECODERS = ("viterbi", "posterior")

# Below this the flip prior is numerically indistinguishable from "never".
_MIN_P = 1e-12


@dataclass
class HMMResult:
    """Posterior and MAP decoding of a two-state branch sequence."""

    posterior: np.ndarray  # (n,) P[branch = 1 | entire trace]
    path: np.ndarray  # (n,) int8 Viterbi branch labels
    log_likelihood: float
    p_flip: float
    diagnostics: dict = field(default_factory=dict)


def gaussian_log_emissions(
    iq: np.ndarray,
    mu_a: np.ndarray,
    mu_b: np.ndarray,
    sigma: float,
) -> np.ndarray:
    """Log-likelihood of each sample under each branch.

    Isotropic complex Gaussian: the real and imaginary parts are independent
    with standard deviation ``sigma`` each, so
    ``log p(z | branch) = -log(2*pi*sigma^2) - |z - mu|^2 / (2*sigma^2)``.

    Returns a ``(2, n)`` array; row 0 is branch A, row 1 is branch B.
    """
    iq = np.asarray(iq)
    if sigma <= 0:
        raise ValueError("sigma must be positive")
    norm = -log(2.0 * np.pi * sigma * sigma)
    inv = 0.5 / (sigma * sigma)
    out = np.empty((2, iq.size), dtype=float)
    out[0] = norm - inv * np.abs(iq - np.asarray(mu_a)) ** 2
    out[1] = norm - inv * np.abs(iq - np.asarray(mu_b)) ** 2
    return out


def _scaled_emissions(log_emit: np.ndarray) -> tuple[list, list, float]:
    """Emission likelihoods rescaled per sample to avoid underflow.

    Each column is divided by its own maximum, so both entries land in
    ``(0, 1]``. The discarded per-sample factors are summed and returned so
    the total log-likelihood can be reconstructed.
    """
    peak = log_emit.max(axis=0)
    rel = np.exp(log_emit - peak)
    return rel[0].tolist(), rel[1].tolist(), float(peak.sum())


def forward_backward(
    log_emit: np.ndarray, p_flip: float
) -> tuple[np.ndarray, float]:
    """Posterior P[branch = 1 | all samples] and the total log-likelihood.

    Symmetric two-state chain: the branch persists with probability
    ``1 - p_flip`` between consecutive samples. The recursion runs in the
    scaled (per-sample normalised) domain, which is exact and avoids both
    underflow and a log/exp call per step.
    """
    p = float(min(max(p_flip, _MIN_P), 0.5))
    stay = 1.0 - p
    e0, e1, log_scale = _scaled_emissions(log_emit)
    n = len(e0)

    # Forward pass. a0, a1 are normalised filtered probabilities.
    fwd0 = [0.0] * n
    fwd1 = [0.0] * n
    a0 = 0.5 * e0[0]
    a1 = 0.5 * e1[0]
    s = a0 + a1
    log_scale += log(s)
    a0 /= s
    a1 /= s
    fwd0[0] = a0
    fwd1[0] = a1
    for k in range(1, n):
        b0 = a0 * stay + a1 * p
        b1 = a0 * p + a1 * stay
        a0 = b0 * e0[k]
        a1 = b1 * e1[k]
        s = a0 + a1
        log_scale += log(s)
        a0 /= s
        a1 /= s
        fwd0[k] = a0
        fwd1[k] = a1

    # Backward pass, folding into the posterior as it goes.
    post = np.empty(n, dtype=float)
    b0 = b1 = 1.0
    post[n - 1] = fwd1[n - 1]
    for k in range(n - 2, -1, -1):
        c0 = b0 * e0[k + 1]
        c1 = b1 * e1[k + 1]
        b0 = c0 * stay + c1 * p
        b1 = c0 * p + c1 * stay
        s = b0 + b1
        b0 /= s
        b1 /= s
        g0 = fwd0[k] * b0
        g1 = fwd1[k] * b1
        post[k] = g1 / (g0 + g1)
    return post, log_scale


def viterbi(log_emit: np.ndarray, p_flip: float) -> np.ndarray:
    """Most likely branch sequence (MAP path) as an int8 array."""
    p = float(min(max(p_flip, _MIN_P), 0.5))
    log_stay = log(1.0 - p)
    log_flip = log(p)
    l0 = log_emit[0].tolist()
    l1 = log_emit[1].tolist()
    n = len(l0)

    back0 = [0] * n
    back1 = [0] * n
    d0 = log(0.5) + l0[0]
    d1 = log(0.5) + l1[0]
    for k in range(1, n):
        stay0 = d0 + log_stay
        from1 = d1 + log_flip
        if stay0 >= from1:
            n0, back0[k] = stay0, 0
        else:
            n0, back0[k] = from1, 1
        stay1 = d1 + log_stay
        from0 = d0 + log_flip
        if stay1 >= from0:
            n1, back1[k] = stay1, 1
        else:
            n1, back1[k] = from0, 0
        d0 = n0 + l0[k]
        d1 = n1 + l1[k]

    path = np.empty(n, dtype=np.int8)
    state = 1 if d1 > d0 else 0
    path[n - 1] = state
    for k in range(n - 1, 0, -1):
        state = back1[k] if state else back0[k]
        path[k - 1] = state
    return path


def decode(
    iq: np.ndarray,
    mu_a: np.ndarray,
    mu_b: np.ndarray,
    sigma: float,
    p_flip: float,
    refine_rate: int = 0,
) -> HMMResult:
    """Full two-state decode: posterior + Viterbi path.

    ``refine_rate`` runs that many Baum-Welch updates of ``p_flip`` alone
    (the emission model is held fixed, since it is learned separately and
    is far better constrained than a single scalar rate).
    """
    log_emit = gaussian_log_emissions(iq, mu_a, mu_b, sigma)
    p = float(p_flip)
    loglik = float("nan")
    history = [p]
    for _ in range(max(0, int(refine_rate))):
        post, loglik = forward_backward(log_emit, p)
        # Expected number of switches, approximated from the marginals. Exact
        # pairwise posteriors would need a second sweep; the marginal estimate
        # is unbiased enough for a rate that only sets the detector's prior.
        p = float(np.mean(np.abs(np.diff(post))))
        p = min(max(p, _MIN_P), 0.5)
        history.append(p)
    post, loglik = forward_backward(log_emit, p)
    path = viterbi(log_emit, p)
    return HMMResult(
        posterior=post,
        path=path,
        log_likelihood=loglik,
        p_flip=p,
        diagnostics={"p_flip_history": history},
    )


def decode_with_rate(
    x: np.ndarray,
    mu_a: np.ndarray,
    mu_b: np.ndarray,
    sigma: float,
    p_flip_init: float = 1e-3,
    n_iter: int = 4,
) -> tuple[HMMResult, float, list]:
    """Decode while re-estimating the flip rate from the data (hard EM).

    The rate only sets the decoder's prior, so a handful of hard-EM passes is
    plenty; a full Baum-Welch update would cost a second sweep per iteration for
    no measurable gain here. The emissions are built once, and the iterations
    run forward-backward only -- the Viterbi path is needed just once, at the
    end, to enumerate the events.

    Returns ``(result, p_flip, history)``.
    """
    n = x.size
    log_emit = gaussian_log_emissions(x, mu_a, mu_b, sigma)
    p = float(p_flip_init)
    history = [p]
    for _ in range(max(1, int(n_iter))):
        post, _ = forward_backward(log_emit, p)
        # Thresholded posterior stands in for the Viterbi path here: it is a
        # transition count, and the two agree except inside low-contrast
        # windows where the posterior sits near 1/2 either way.
        crossings = np.count_nonzero(np.diff(post > 0.5))
        p_new = max(crossings / max(n - 1, 1), 1e-9)
        history.append(p_new)
        converged = abs(p_new - p) < 1e-3 * p
        p = p_new
        if converged:
            break
    post, loglik = forward_backward(log_emit, p)
    path = viterbi(log_emit, p)
    res = HMMResult(posterior=post, path=path, log_likelihood=loglik, p_flip=p,
                    diagnostics={"p_flip_history": history})
    return res, p, history


def decoded_path(result: HMMResult, decoder: str = "viterbi") -> np.ndarray:
    """Branch sequence to read events off, under one of two decoding rules.

    Both rules use the *same* fitted emissions and the same flip prior --
    :func:`decode_with_rate` estimates the rate from the posterior either way
    -- so choosing between them isolates the decoding step and nothing else.

    ``"viterbi"``
        The single most likely branch *sequence* (MAP path). It is a jointly
        valid trajectory, and the transition prior is paid once per flip, so a
        brief dip in the evidence does not create a flip *pair* unless it earns
        its keep globally. Conservative: fewer, cleaner transitions.

    ``"posterior"``
        The per-sample marginal thresholded at 1/2 (MAP marginal, from
        forward-backward). It maximises the expected number of correctly
        labelled *samples*, which is a different objective -- and the resulting
        sequence need not be a likely trajectory at all, since nothing couples
        neighbouring decisions once the marginals are formed. It will follow a
        short excursion that Viterbi smooths over, so it tends to recover more
        real flips and also to invent more.

    Which is better is an empirical question about the operating point, not a
    matter of one dominating the other; see ``docs/reconstruction.md`` §12e and
    ``notebooks/reconstruction_evaluation.ipynb``.
    """
    if decoder == "viterbi":
        return np.asarray(result.path)
    if decoder == "posterior":
        return (np.asarray(result.posterior, dtype=float) > 0.5).astype(np.int8)
    raise ValueError(
        f"decoder must be one of {DECODERS}; got {decoder!r}")


def forward_backward_reference(
    log_emit: np.ndarray, log_transition: np.ndarray
) -> tuple[np.ndarray, float]:
    """Transparent general-K log-domain forward-backward, for testing.

    ``log_transition[i, j]`` is log P(state j at k+1 | state i at k). Returns
    the ``(K, n)`` posterior and the total log-likelihood. Deliberately
    written for clarity, not speed.
    """
    from scipy.special import logsumexp

    k_states, n = log_emit.shape
    alpha = np.full((k_states, n), -np.inf)
    alpha[:, 0] = -log(k_states) + log_emit[:, 0]
    for k in range(1, n):
        alpha[:, k] = (
            logsumexp(alpha[:, k - 1][:, None] + log_transition, axis=0)
            + log_emit[:, k]
        )
    beta = np.zeros((k_states, n))
    for k in range(n - 2, -1, -1):
        beta[:, k] = logsumexp(
            log_transition + (log_emit[:, k + 1] + beta[:, k + 1])[None, :],
            axis=1,
        )
    loglik = float(logsumexp(alpha[:, -1]))
    log_gamma = alpha + beta - loglik
    return np.exp(log_gamma), loglik
