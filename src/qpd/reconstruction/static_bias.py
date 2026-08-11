"""Parity-flip reconstruction at a fixed offset charge.

When ``n_g`` is held constant the readout collapses to the textbook case: the
two parity branches are two *stationary* points in the I/Q plane, and the trace
is a random telegraph hopping between them under additive Gaussian noise. None
of the machinery the ramped case needs applies -- there is no fold period, no
parity-blind point sweeping past, no sawtooth reset swapping the branch labels,
and no ramp phase for a charge jump to shift. So this is deliberately a separate
and much shorter entry point rather than a special case of
:func:`~qpd.reconstruction.reconstruct.reconstruct_parity_flips_ramped`.

What is left is a two-component Gaussian mixture in the plane plus a dwell-time
prior. The mixture is fitted blind: the two blob centres differ along a single
line, so the trace's principal axis *is* the discrimination axis, and projecting
onto it turns the problem into 1-D. The projected histogram is bimodal, and its
two modes and common width are recovered by a short EM. The branch sequence then
comes from the same two-state HMM used in the ramped case.

Two things are worth knowing about this regime:

* **The contrast is whatever the bias point gives you, and it does not vary.**
  In the ramped case every flip eventually lands somewhere with good separation;
  here, a bias point near the parity-blind charge (n_g = 0.25 mod 0.5) yields
  two blobs that overlap and *no* amount of decoding recovers the flips. The
  fitted separation is reported so this is visible rather than silent.
* **Charge jumps are not nuisances here -- they are the only thing that moves
  the blobs.** A jump shifts both centres to a new pair of positions and the
  fitted mixture goes stale. ``segment_blocks`` re-fits the mixture in blocks so
  a drifting or jumping bias does not quietly poison the whole trace; with a
  genuinely constant bias the default single block is correct and cheapest.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .emission import validate_trace
from .events import extract_flips
from .hmm import decode_with_rate, decoded_path

__all__ = ["StaticBlobModel", "StaticReconstructionResult", "fit_two_blobs",
           "reconstruct_parity_flips_static"]


@dataclass
class StaticBlobModel:
    """Two stationary branch means and a shared noise width, on one axis."""

    direction: complex  # unit vector of the discrimination axis in I/Q
    origin: complex  # I/Q point the projection is measured from
    mu_a: float  # projected centre of branch A
    mu_b: float  # projected centre of branch B
    sigma: float  # per-quadrature noise std (shared by both blobs)
    weight_a: float  # mixture weight of branch A (steady-state occupancy)
    n_iter: int = 0
    diagnostics: dict = field(default_factory=dict)

    def project(self, iq: np.ndarray) -> np.ndarray:
        """Map complex I/Q onto the real discrimination axis."""
        return np.real(np.conj(self.direction) * (np.asarray(iq) - self.origin))

    @property
    def separation(self) -> float:
        """Distance between the two blob centres, on the projected axis."""
        return abs(self.mu_a - self.mu_b)

    @property
    def contrast(self) -> float:
        """Blob separation in units of the noise std.

        The single number that decides whether this bias point is usable. Below
        about 1 the blobs overlap badly and flips are not recoverable; the
        ramped case has no equivalent because its contrast sweeps.
        """
        return self.separation / self.sigma if self.sigma > 0 else np.inf


def fit_two_blobs(
    iq: np.ndarray,
    n_iter: int = 60,
    tol: float = 1e-10,
) -> StaticBlobModel:
    """Fit two equal-width Gaussian blobs to a static-bias I/Q trace, blind.

    The two branch means differ along one line, so the cloud's principal axis is
    the discrimination axis and the projection onto it carries all the parity
    information. On that axis the data is a two-component Gaussian mixture with a
    shared width, fitted here by EM.

    EM is initialised by splitting the projection at its median, which is a
    deliberate choice: a k-means-style split at the extremes biases the centres
    outward when the blobs overlap, whereas the median split is unbiased for the
    balanced occupancy this process has (equal dwell times in the two states).
    """
    iq = np.asarray(iq)
    if iq.size < 8:
        raise ValueError("need at least 8 samples to fit two blobs")

    origin = complex(np.mean(iq))
    d = iq - origin
    x_, y_ = d.real, d.imag
    cov = np.array([[np.mean(x_ * x_), np.mean(x_ * y_)],
                    [np.mean(x_ * y_), np.mean(y_ * y_)]])
    vals, vecs = np.linalg.eigh(cov)  # ascending
    major = vecs[:, 1]
    direction = complex(major[0], major[1])
    # The minor axis carries noise only, so it is an independent read of sigma
    # that does not depend on the mixture fit converging.
    sigma_minor = float(np.sqrt(max(vals[0], 1e-30)))

    x = np.real(np.conj(direction) * d)
    split = float(np.median(x))
    mu_a = float(np.mean(x[x >= split])) if np.any(x >= split) else split
    mu_b = float(np.mean(x[x < split])) if np.any(x < split) else split
    sigma = max(sigma_minor, 1e-12)
    w_a = 0.5

    loglik_prev = -np.inf
    it = 0
    for it in range(1, int(n_iter) + 1):
        inv = 0.5 / (sigma * sigma)
        la = np.log(max(w_a, 1e-300)) - inv * (x - mu_a) ** 2
        lb = np.log(max(1.0 - w_a, 1e-300)) - inv * (x - mu_b) ** 2
        m = np.maximum(la, lb)
        tot = m + np.log(np.exp(la - m) + np.exp(lb - m))
        resp = np.exp(la - tot)  # P[branch A | sample]

        n_a = float(resp.sum())
        n_b = float(x.size - n_a)
        if n_a < 1e-6 or n_b < 1e-6:
            break
        mu_a = float((resp * x).sum() / n_a)
        mu_b = float(((1.0 - resp) * x).sum() / n_b)
        var = float((resp * (x - mu_a) ** 2
                     + (1.0 - resp) * (x - mu_b) ** 2).sum() / x.size)
        sigma = float(np.sqrt(max(var, 1e-30)))
        w_a = n_a / x.size

        loglik = float(tot.sum() - 0.5 * x.size * np.log(2 * np.pi * sigma ** 2))
        if abs(loglik - loglik_prev) <= tol * max(abs(loglik), 1.0):
            break
        loglik_prev = loglik

    # Single-Gaussian reference, for the degeneracy test below. Fitting two
    # components to one blob always "succeeds" -- EM will split a single
    # Gaussian into a spurious pair -- so the fitted separation on its own
    # cannot say whether the bias point carries any parity information.
    sigma_one = float(np.std(x))
    loglik_one = float(-0.5 * x.size * (1.0 + np.log(2 * np.pi * sigma_one ** 2)))
    loglik_two = loglik_prev if np.isfinite(loglik_prev) else loglik_one

    # Order the blobs so A is the upper one on the axis; the labelling is
    # arbitrary (a global even/odd relabel) and does not affect flip timing.
    if mu_b > mu_a:
        mu_a, mu_b = mu_b, mu_a
        w_a = 1.0 - w_a

    return StaticBlobModel(
        direction=direction, origin=origin, mu_a=mu_a, mu_b=mu_b,
        sigma=sigma, weight_a=w_a, n_iter=it,
        diagnostics={
            "sigma_minor_axis": sigma_minor,
            "sigma_ratio_fit_over_minor": sigma / max(sigma_minor, 1e-30),
            "separation": abs(mu_a - mu_b),
            "contrast": abs(mu_a - mu_b) / max(sigma, 1e-30),
            "em_iterations": it,
            "loglik_two_blob": loglik_two,
            "loglik_one_blob": loglik_one,
            "loglik_gain": loglik_two - loglik_one,
        },
    )


@dataclass
class StaticReconstructionResult:
    """Reconstructed flip times at a fixed bias, and the model behind them."""

    flip_times: np.ndarray  # [s] estimated flip times
    confidence: np.ndarray  # [0, 1] posterior swing at each flip
    posterior: np.ndarray  # (n,) P[branch B | trace]
    branch: np.ndarray  # (n,) int8 Viterbi branch labels
    model: StaticBlobModel | None = None
    p_flip: float = 0.0
    diagnostics: dict = field(default_factory=dict)

    @property
    def rate_hz(self) -> float:
        """Estimated per-state tunnelling rate."""
        return self.diagnostics.get("rate_hz", float("nan"))

    @property
    def contrast(self) -> float:
        """Blob separation in units of the noise std."""
        return self.model.contrast if self.model is not None else float("nan")

    @property
    def sample_fidelity(self) -> float:
        """Single-sample parity assignment fidelity, from the fitted contrast.

        This is the conventional readout-fidelity figure: with two equal-width
        Gaussians a contrast ``C`` apart and the threshold at their midpoint,
        each sample is misassigned with probability
        ``0.5 * erfc(C / (2*sqrt(2)))``, and the fidelity is one minus that.

        It describes *one sample in isolation*, so it is a pessimistic view of
        what this pipeline achieves: the decoder does not classify samples
        independently, it integrates over a dwell (see
        :attr:`decoded_fidelity`). Quote this when comparing against a
        single-shot readout, not when describing the reconstruction.
        """
        from scipy.special import erfc
        c = self.contrast
        if not np.isfinite(c):
            return float("nan")
        return float(1.0 - 0.5 * erfc(c / (2.0 * np.sqrt(2.0))))

    @property
    def decoded_fidelity(self) -> float:
        """Expected fraction of samples assigned to the correct branch.

        Read straight off the posterior: a sample whose posterior is ``g`` is
        misassigned with probability ``min(g, 1-g)`` *under the fitted model*,
        so the mean of that over the trace is the model's own estimate of its
        per-sample error rate. Because the HMM pools evidence across a whole
        dwell, this is far higher than :attr:`sample_fidelity`.

        **It is conditional on the model being right.** At a parity-blind bias
        the fitted two-blob model is spurious, and the decoder can be
        confidently wrong -- this number stays high while the output is
        meaningless. Always read it together with :attr:`degenerate` and
        :attr:`contrast`; it measures self-consistency, not correctness.
        """
        g = np.asarray(self.posterior, dtype=float)
        if g.size == 0:
            return float("nan")
        return float(1.0 - np.mean(np.minimum(g, 1.0 - g)))

    @property
    def degenerate(self) -> bool:
        """True when the bias point carries no usable parity information.

        Near the parity-blind charge (n_g = 0.25 mod 0.5) the two branches
        coincide, and the fitted separation is *not* a reliable warning on its
        own: EM splits a single noise blob into a spurious pair, reporting a
        contrast near 1 rather than 0. The flag therefore requires *both* a low
        integrated detectability (contrast x sqrt(dwell), which collapses when
        the decoder is segmenting noise) and a low per-sample contrast -- the
        second condition being what stops a genuinely fast telegraph, whose
        dwell is legitimately short, from being mislabelled.

        When this is True the ``flip_times`` are meaningless -- move the bias
        point rather than trying to salvage them.
        """
        return bool(self.diagnostics.get("degenerate", False))


def reconstruct_parity_flips_static(
    iq: np.ndarray,
    sample_rate: float,
    t0: float = 0.0,
    p_flip_init: float = 1e-3,
    n_rate_iterations: int = 4,
    min_confidence: float = 0.0,
    segment_blocks: int = 1,
    min_detectability: float = 70.0,
    min_contrast: float = 1.5,
    model: StaticBlobModel | None = None,
    decoder: str = "viterbi",
) -> StaticReconstructionResult:
    """Recover parity-flip times from a **fixed-offset-charge** trace, blind.

    Use this when ``n_g`` is held constant, so the readout is two stationary
    blobs in the I/Q plane with telegraph switching between them. For a swept
    ``n_g`` use
    :func:`~qpd.reconstruction.reconstruct.reconstruct_parity_flips_ramped`
    instead -- the moving branch means, the parity-blind crossings, and the
    sawtooth resets it handles are all absent here, and this routine would model
    them as noise.

    Pipeline: fit two equal-width Gaussian blobs on the trace's principal axis
    (:func:`fit_two_blobs`) -> two-state HMM decode with the flip rate estimated
    from the data -> read flips off the transitions.

    Parameters
    ----------
    iq : complex array
        Readout trace, ``I + 1j*Q``, on a uniform time grid.
    sample_rate : float
        Sampling rate [Hz].
    t0 : float
        Time of the first sample [s]; flip times are reported on this axis.
    p_flip_init : float
        Starting per-sample flip probability; re-estimated from the data, so it
        only has to be within an order of magnitude.
    n_rate_iterations : int
        Hard-EM refinements of the flip rate.
    min_confidence : float
        Drop flips whose posterior swing is below this. Raising it trades recall
        for precision, which matters most at low contrast.
    segment_blocks : int
        Fit the blob model independently in this many equal blocks. The default
        of 1 is right for a genuinely constant bias. Raise it if the bias point
        drifts or takes a charge jump mid-trace: unlike the ramped case, a jump
        here moves the blobs outright, and one stale global fit would corrupt the
        whole decode.
    min_detectability : float
        Threshold on ``contrast * sqrt(dwell in samples)`` -- the contrast
        integrated over one dwell. Only declares the bias point unusable when
        ``min_contrast`` also fails; see
        :attr:`StaticReconstructionResult.degenerate`.
    min_contrast : float
        Per-sample contrast above which the bias point is accepted regardless of
        dwell. This is what keeps a genuinely *fast* telegraph from being
        mislabelled: a device tunnelling at 3 kHz has a short dwell and so a low
        integrated detectability, but with the branches well separated per
        sample its flips are still recovered (measured F1 0.97 at 3 kHz on the
        reference device).
    model : StaticBlobModel, optional
        Pre-fitted blob model; fitted from ``iq`` when omitted. Supplying one is
        the way to force a known bias point, or to reuse a fit across traces.
    decoder : {"viterbi", "posterior"}
        Which decoding rule turns the HMM output into events -- the most likely
        branch *sequence*, or the per-sample marginal thresholded at 1/2. Both
        run on the same emissions and the same flip prior (the rate is
        estimated from the posterior either way), so this switches the decoding
        step alone. See :func:`~qpd.reconstruction.hmm.decoded_path`, and
        ``notebooks/reconstruction_evaluation.ipynb`` for the measured
        comparison.

    Returns
    -------
    StaticReconstructionResult
    """
    iq = validate_trace(iq, sample_rate, min_samples=8)
    dt = 1.0 / float(sample_rate)
    n = iq.size
    blocks = max(1, int(segment_blocks))

    if model is not None or blocks == 1:
        fitted = model if model is not None else fit_two_blobs(iq)
        x = fitted.project(iq)
        mu_a = np.full(n, fitted.mu_a)
        mu_b = np.full(n, fitted.mu_b)
        sigma = fitted.sigma
        per_block = [fitted]
    else:
        # Per-block fits, evaluated on a common axis so the decoded branch label
        # stays comparable across block boundaries.
        fitted = fit_two_blobs(iq)
        x = fitted.project(iq)
        mu_a = np.empty(n)
        mu_b = np.empty(n)
        sigmas = []
        per_block = []
        edges = np.linspace(0, n, blocks + 1).astype(int)
        for lo, hi in zip(edges[:-1], edges[1:]):
            if hi - lo < 8:
                mu_a[lo:hi], mu_b[lo:hi] = fitted.mu_a, fitted.mu_b
                sigmas.append(fitted.sigma)
                per_block.append(fitted)
                continue
            b = fit_two_blobs(iq[lo:hi])
            # A block holding too few flips lets EM converge on a degenerate
            # split, whose centres are nowhere near the true branches. Falling
            # back to the global fit for that block is always safe: the bias is
            # constant by assumption, so the global centres are valid
            # everywhere, and only the *local* refinement is lost.
            if b.contrast < 0.5 * fitted.contrast:
                mu_a[lo:hi], mu_b[lo:hi] = fitted.mu_a, fitted.mu_b
                sigmas.append(fitted.sigma)
                per_block.append(fitted)
                continue
            xb = b.project(iq[lo:hi])
            # Re-express the block's centres on the global axis by matching the
            # block's own projected means to their global-axis counterparts.
            hi_mask = xb >= 0.5 * (b.mu_a + b.mu_b)
            seg = x[lo:hi]
            m_hi = float(np.mean(seg[hi_mask])) if np.any(hi_mask) else fitted.mu_a
            m_lo = float(np.mean(seg[~hi_mask])) if np.any(~hi_mask) else fitted.mu_b
            mu_a[lo:hi], mu_b[lo:hi] = m_hi, m_lo
            sigmas.append(b.sigma)
            per_block.append(b)
        sigma = float(np.median(sigmas))

    res, p, history = decode_with_rate(x, mu_a, mu_b, sigma, p_flip_init,
                                       n_rate_iterations)
    path = decoded_path(res, decoder)
    times, conf = extract_flips(path, res.posterior, dt, t0=t0)
    if min_confidence > 0:
        keep = conf >= min_confidence
        times, conf = times[keep], conf[keep]

    # Self-consistency test for the bias point. Two statistics are needed
    # because each fails on a different axis, and the failures do not overlap.
    #
    # `detectability` = contrast * sqrt(dwell) separates a parity-blind bias
    # from a usable one at a fixed flip rate: decoding noise inflates the rate
    # and collapses the dwell, so the product falls away sharply (measured on
    # the reference device: 105 at n_g = 0.22, which works, against 6 at
    # n_g = 0.24, which does not). But it is not scale-free -- a genuinely fast
    # telegraph also has a short dwell, and a 3 kHz device scores only ~22 while
    # still reconstructing at F1 0.97.
    #
    # The per-sample `contrast` covers exactly that gap, and fails where
    # detectability works: EM splits a single blob into a spurious pair with
    # contrast ~0.9, so on its own it cannot tell a blind point (0.90) from a
    # marginal but usable one (0.96).
    #
    # Requiring *both* to fail is what makes the flag trustworthy in each
    # direction. The likelihood gain of two blobs over one was also tried and
    # is useless here: it is ~0 for the usable n_g = 0.22 and the blind
    # n_g = 0.24 alike.
    dwell_samples = 1.0 / max(p, 1e-12)
    detectability = fitted.contrast * np.sqrt(dwell_samples)
    degenerate = bool(detectability < float(min_detectability)
                      and fitted.contrast < float(min_contrast))

    return StaticReconstructionResult(
        flip_times=times,
        confidence=conf,
        posterior=res.posterior,
        branch=path,
        model=fitted,
        p_flip=p,
        diagnostics={
            "decoder": decoder,
            "viterbi_path": res.path,
            "rate_hz": p / dt,
            "p_flip_history": history,
            "n_flips": int(times.size),
            "log_likelihood": res.log_likelihood,
            "contrast": fitted.contrast,
            "separation": fitted.separation,
            "sigma": sigma,
            "weight_a": fitted.weight_a,
            "segment_blocks": blocks,
            "degenerate": degenerate,
            "detectability": float(detectability),
            "dwell_samples": float(dwell_samples),
            "loglik_gain_two_vs_one_blob": fitted.diagnostics.get("loglik_gain"),
            "block_models": per_block if blocks > 1 else None,
        },
    )
