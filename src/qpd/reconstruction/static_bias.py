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

from .burst_hmm import decode_burst_aware
from .emission import validate_trace
from .events import extract_flips
from .hmm import decode_with_rate, decoded_path

__all__ = ["StaticBlobModel", "StaticReconstructionResult", "fit_two_blobs",
           "reconstruct_parity_flips_static"]

# More accepted charge events than this in one trace is not physics: real
# jump rates are ~0.1 Hz, so even a 5 s trace should carry 0 or 1 (a few for
# a correlated impact cluster). Beyond it the stationary-blob model itself is
# in doubt -- drift, transients, gain wander -- and acting on the boundaries
# would shred good data, so segmentation is abandoned (and said so) instead.
_MAX_PLAUSIBLE_EVENTS = 3

# Dead-time test for post-jump segments. The global `min_detectability` (70)
# is calibrated for whole-trace flagging and is far stricter than what
# separates a *disastrous* segment from a marginal one: a genuinely blind
# landing makes the decoder toggle, which collapses the decoded dwell, and
# with it contrast * sqrt(dwell) -- measured 6.8 at the parity-blind landing
# against 15.8 at the nearest marginal-but-decodable bias (10 kSa/s, 17 Hz).
# 12 splits the two with margin both ways. A marginal segment is *decoded and
# reported* (same semantics as the global degenerate flag, which flags but
# still returns); only the toggle-collapse case becomes dead time.
_SEGMENT_DEAD_DETECTABILITY = 12.0


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
    # Detected offset-charge events [s], same axis as flip_times. Mirrors the
    # ramped ReconstructionResult attribute of the same name.
    charge_jump_times: np.ndarray = field(
        default_factory=lambda: np.empty(0))

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
    burst_aware: bool = True,
    burst_rate_hz: float = 1.0,
    p_burst: float = 0.3,
    burst_tau: float = 1e-3,
    segment_charge_jumps: bool = False,
    charge_min_gain_nats: float = 15.0,
    charge_min_segment_samples: int = 2500,
    boundary_guard_samples: int = 3,
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
        of 1 is right for a genuinely constant bias; raise it for slow *drift*.
        Discrete charge jumps are handled by ``segment_charge_jumps`` (the
        detected boundaries supersede this grid when both are active, since
        the measured edges answer the same question better than equal-width
        blocks).
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
    burst_aware : bool
        Decode with the parity x regime chain of
        :mod:`~qpd.reconstruction.burst_hmm` (the default) instead of the
        plain two-state HMM. The global flip prior smooths quasiparticle
        bursts down to two or three recovered flips regardless of their true
        size (issue #40); the burst-aware chain gives burst windows their own
        flip probability, at no change to background behaviour, which is why
        it is on by default. The reported ``rate_hz`` then refers to the
        *quiet* regime, i.e. it is the background rate no longer inflated by
        the bursts. Pass ``False`` for the two-state decoder -- it is ~6x
        cheaper and identical on burst-free traces.
    burst_rate_hz : float
        Expected burst occurrence rate [Hz], only used when ``burst_aware``.
        Sets the regime entry prior; logarithmic sensitivity, so an order of
        magnitude either way moves the burst detection threshold by about one
        flip of evidence.
    p_burst : float
        Pinned per-sample flip probability inside a burst (not fitted).
    burst_tau : float
        Burst decay time [s]; sets the regime exit prior.
    segment_charge_jumps : bool
        Detect offset-charge events and refit the blob model per segment.
        Unlike the ramped case, a jump here moves the blobs outright: an
        undetected one smears flips, silently kills the post-jump stretch
        when the new bias lands near the parity-blind charge (measured F1
        0.95 -> 0.30 with the reported rate inflated 16 -> 40 Hz), or
        starves a coincident burst -- and real impact events produce charge
        jumps *and* bursts together, so that failure lands on the signal.
        Detection is blind
        (:func:`~qpd.reconstruction.charge_events.detect_charge_events`);
        on a trace with no verified event the decode is identical to
        ``segment_charge_jumps=False``, and a post-jump segment whose refit
        fails the toggle-collapse test is declared **dead time** --
        reported, masked, never silently decoded with a stale model.

        **Off by default, deliberately.** On measured devices the emission
        model is not piecewise-constant: 1/f charge noise is physically a
        dense train of micro-jumps, and pulsed calibration sources (LED
        flashes) are genuine sharp model steps -- on the reference LED
        dataset the detector locks onto the 60 ms flash comb with formally
        enormous significance. Both are correct detections of the wrong
        thing, and no within-trace statistic can separate them from the
        rare discrete jump this option targets (several were tried and
        measured; see ``docs/reconstruction.md`` section 12g). Turn it on
        for dark data at a stable bias, or when a specific trace is
        suspected; read ``diagnostics["charge_event_times"]`` and
        ``live_fraction`` before trusting the output either way.
    charge_min_gain_nats : float
        Acceptance threshold on the detector's verified likelihood ratio.
    charge_min_segment_samples : int
        Shortest segment that gets its own refit; anything shorter between
        detected boundaries (or against the trace ends) is dead time.
    boundary_guard_samples : int
        Drop flips within this many samples of a detected boundary: branch
        identity does not survive a jump, so a transition there is
        bookkeeping, not physics. (A quasiparticle-tunnelling-induced jump
        does flip the parity, but the flip cannot be counted honestly.)

    Returns
    -------
    StaticReconstructionResult
    """
    iq = validate_trace(iq, sample_rate, min_samples=8)
    dt = 1.0 / float(sample_rate)
    n = iq.size
    blocks = max(1, int(segment_blocks))

    fitted = model if model is not None else fit_two_blobs(iq)
    x = fitted.project(iq)

    # Charge-event detection runs first, on the global fit -- also when a
    # model was supplied, since a forced bias point can still jump.
    events = None
    if segment_charge_jumps and n >= 2 * int(charge_min_segment_samples):
        # Local import: charge_events imports StaticBlobModel from here.
        from .charge_events import detect_charge_events
        events = detect_charge_events(
            iq, fitted, dt,
            min_gain_nats=float(charge_min_gain_nats),
            min_segment_samples=int(charge_min_segment_samples),
            p_flip=p_flip_init)

    abandoned = None
    if events is not None and events.boundaries.size > _MAX_PLAUSIBLE_EVENTS:
        abandoned, events = events, None

    live = None
    seg_edges = None
    seg_models: list[StaticBlobModel | None] | None = None
    dead_idx: list[tuple[int, int]] = []
    if events is not None and events.boundaries.size:
        # Verified boundaries supersede the equal-width segment_blocks grid:
        # both mechanisms answer "where does the blob model change", and the
        # detected edges are the measured answer.
        blocks = 1
        min_seg = int(charge_min_segment_samples)
        seg_edges = np.concatenate(([0], events.boundaries, [n])).astype(int)
        mu_a = np.empty(n)
        mu_b = np.empty(n)
        sigmas = []
        seg_models = []
        live = np.ones(n, dtype=bool)
        per_block = [fitted]
        for lo, hi in zip(seg_edges[:-1], seg_edges[1:]):
            seg = x[lo:hi]
            seg_fit = None
            m_hi = m_lo = 0.0
            dead = hi - lo < min_seg
            if not dead:
                seg_fit = fit_two_blobs(iq[lo:hi])
                # Re-express the segment's centres on the global axis
                # *geometrically*: map each fitted centre back to the I/Q
                # plane and project it. The segment_blocks path below instead
                # averages hard-classified samples, whose conditional means
                # are truncation-biased apart at contrast ~1 -- enough to
                # make a marginal post-jump segment decode as a toggling mess
                # and fail the dead test it should pass.
                for_c = np.conj(fitted.direction)
                pa, pb = (float(np.real(for_c * (
                    seg_fit.origin + seg_fit.direction * mu - fitted.origin)))
                    for mu in (seg_fit.mu_a, seg_fit.mu_b))
                m_hi, m_lo = (pa, pb) if pa >= pb else (pb, pa)
                # Per-segment health: the same *shape* of rule as the global
                # `degenerate` flag but with the segment-calibrated
                # detectability cut (see _SEGMENT_DEAD_DETECTABILITY),
                # evaluated in the segment's own frame. There is deliberately
                # NO fallback to the global centres here -- after a jump they
                # are exactly the wrong answer, and installing them would
                # reintroduce the silent failure this whole branch exists to
                # close.
                xb = seg_fit.project(iq[lo:hi])
                _, p_seg, _ = decode_with_rate(
                    xb, np.full(hi - lo, seg_fit.mu_a),
                    np.full(hi - lo, seg_fit.mu_b),
                    seg_fit.sigma, p_flip_init, 2)
                det_seg = seg_fit.contrast * np.sqrt(1.0 / max(p_seg, 1e-12))
                dead = bool(det_seg < _SEGMENT_DEAD_DETECTABILITY
                            and seg_fit.contrast < float(min_contrast))
            if dead:
                # Equalized means carry no parity evidence: the decoder
                # coasts through on the transition prior alone, which is the
                # honest treatment of a stretch whose emissions are unusable.
                mu_a[lo:hi] = mu_b[lo:hi] = float(np.mean(seg))
                live[lo:hi] = False
                dead_idx.append((int(lo), int(hi)))
                seg_models.append(None)
            else:
                mu_a[lo:hi], mu_b[lo:hi] = m_hi, m_lo
                sigmas.append(seg_fit.sigma)
                seg_models.append(seg_fit)
        sigma = float(np.median(sigmas)) if sigmas else fitted.sigma
        if bool(live.all()):
            live = None
    elif model is not None or blocks == 1:
        mu_a = np.full(n, fitted.mu_a)
        mu_b = np.full(n, fitted.mu_b)
        sigma = fitted.sigma
        per_block = [fitted]
    else:
        # Per-block fits, evaluated on a common axis so the decoded branch label
        # stays comparable across block boundaries.
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

    if burst_aware:
        bres = decode_burst_aware(
            x, mu_a, mu_b, sigma, dt, p_flip_init=p_flip_init,
            n_iter=n_rate_iterations, p_burst=p_burst,
            burst_rate_hz=burst_rate_hz, burst_tau=burst_tau, live=live)
        p = bres.p_quiet
        history = bres.diagnostics["p_quiet_history"]
        posterior = bres.parity_posterior
        viterbi_path = bres.parity_path
        if decoder == "viterbi":
            path = viterbi_path
        elif decoder == "posterior":
            path = (posterior > 0.5).astype(np.int8)
        else:
            raise ValueError(
                f"decoder must be 'viterbi' or 'posterior'; got {decoder!r}")
        loglik = bres.log_likelihood
    else:
        res, p, history = decode_with_rate(x, mu_a, mu_b, sigma, p_flip_init,
                                           n_rate_iterations)
        path = decoded_path(res, decoder)
        posterior = res.posterior
        viterbi_path = res.path
        loglik = res.log_likelihood
        bres = None
    times, conf = extract_flips(path, posterior, dt, t0=t0)
    n_suppressed = 0
    if events is not None and events.boundaries.size:
        # Flips at a boundary are bookkeeping (branch identity does not
        # survive a jump), and flips inside dead windows are decoded against
        # equalized means, i.e. noise. Both go, before the confidence cut so
        # the reported counts compose.
        keep = np.ones(times.size, dtype=bool)
        guard = float(boundary_guard_samples) * dt
        for b in events.boundaries:
            keep &= np.abs(times - (t0 + b * dt)) > guard
        for lo, hi in dead_idx:
            keep &= ~((times >= t0 + lo * dt) & (times < t0 + hi * dt))
        n_suppressed = int(np.count_nonzero(~keep))
        times, conf = times[keep], conf[keep]
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
    # When the trace was segmented, the global fit is a blend across the jump
    # and its contrast understates the usable segments; judge health on the
    # longest live segment instead. No live segment at all is degenerate by
    # definition.
    health = fitted
    no_live = False
    if seg_models is not None:
        spans = [(hi - lo, m) for (lo, hi), m in
                 zip(zip(seg_edges[:-1], seg_edges[1:]), seg_models)
                 if m is not None]
        if spans:
            health = max(spans, key=lambda s: s[0])[1]
        else:
            no_live = True
    dwell_samples = 1.0 / max(p, 1e-12)
    detectability = health.contrast * np.sqrt(dwell_samples)
    degenerate = bool(no_live
                      or (detectability < float(min_detectability)
                          and health.contrast < float(min_contrast)))

    charge_extra = {}
    jump_times = np.empty(0)
    if abandoned is not None:
        charge_extra = {
            "segment_charge_jumps": True,
            "charge_segmentation_abandoned": True,
            "charge_event_times": t0 + abandoned.boundaries * dt,
            "charge_event_gain_nats": abandoned.gain_nats,
            "charge_scan": abandoned.scan,
        }
    if events is not None:
        jump_times = t0 + events.boundaries * dt
        charge_extra = {
            "segment_charge_jumps": True,
            "charge_event_times": jump_times,
            "charge_event_gain_nats": events.gain_nats,
            "charge_event_localization_s": events.localization_sigma * dt,
            "charge_scan": events.scan,
        }
        if events.boundaries.size:
            charge_extra.update({
                "segment_edges": seg_edges,
                "segment_models": seg_models,
                "dead_windows": [(t0 + lo * dt, t0 + hi * dt)
                                 for lo, hi in dead_idx],
                "live_fraction": (float(np.mean(live)) if live is not None
                                  else 1.0),
                "n_boundary_suppressed": n_suppressed,
                "segment_blocks_superseded": int(segment_blocks) > 1,
                "global_model": fitted,
            })

    extra = {}
    if bres is not None:
        # Contiguous burst-regime stretches of the Viterbi path, as time
        # windows on the same axis as the flip times.
        r = np.asarray(bres.regime_path)
        edges = np.flatnonzero(np.diff(r) != 0) + 1
        starts = np.r_[0, edges][np.r_[r[0], r[edges]] == 1]
        ends = np.r_[edges, r.size][np.r_[r[0], r[edges]] == 1]
        extra = {
            "burst_aware": True,
            "burst_posterior": bres.burst_posterior,
            "regime_path": r,
            "burst_windows": [(t0 + lo * dt, t0 + hi * dt)
                              for lo, hi in zip(starts, ends)],
            "p_quiet": bres.p_quiet,
            "p_burst": bres.p_burst,
            "epsilon": bres.diagnostics["epsilon"],
            "p_global_seed": bres.diagnostics["p_global_seed"],
        }

    return StaticReconstructionResult(
        flip_times=times,
        confidence=conf,
        posterior=posterior,
        branch=path,
        # On a segmented trace the global fit is a blend; the longest live
        # segment's model is the representative one (the global fit is kept
        # under diagnostics["global_model"]).
        model=health,
        p_flip=p,
        charge_jump_times=jump_times,
        diagnostics={
            "decoder": decoder,
            "viterbi_path": viterbi_path,
            "rate_hz": p / dt,
            "p_flip_history": history,
            "n_flips": int(times.size),
            "log_likelihood": loglik,
            **extra,
            **charge_extra,
            "contrast": health.contrast,
            "separation": health.separation,
            "sigma": sigma,
            "weight_a": health.weight_a,
            "segment_blocks": blocks,
            "degenerate": degenerate,
            "detectability": float(detectability),
            "dwell_samples": float(dwell_samples),
            "loglik_gain_two_vs_one_blob": health.diagnostics.get("loglik_gain"),
            "block_models": per_block if blocks > 1 else None,
        },
    )
