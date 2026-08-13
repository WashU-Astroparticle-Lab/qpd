"""Calibrate the reconstruction against a *measured* trace, by surrogate replay.

A real trace comes with no truth: nobody knows when the quasiparticles actually
tunnelled, so :func:`~qpd.reconstruction.score_flips` cannot be run on it and
the reconstruction's efficiency and timing accuracy on that specific dataset are
unknown. Every diagnostic the pipeline reports about a real trace
(``contrast``, ``decoded_fidelity``, ``degenerate``) measures *self-consistency
with the fitted model*, not correctness -- and at a bad bias point they stay
healthy while the output is meaningless.

This module closes that gap the only way it can be closed: by simulating traces
that carry the measured trace's own fidelity and *do* have truth.

The chain is

1. **Characterise.** Run the blind reconstruction on the measured trace. What
   comes back is not only a flip list but the fitted forward model -- the
   discrimination axis, the noise sigma, the branch splitting (as a constant at
   fixed bias, or as a profile over ramp phase when ``n_g`` is swept), the ramp
   reset schedule, and the decoded tunnelling rate. That model *is* the
   measurement's fidelity, expressed in the only units that matter to the
   decoder: branch separation in units of noise (:class:`TraceFidelity`).

2. **Replay.** Draw a fresh telegraph trajectory at the measured rate and push
   it through the fitted model with fresh noise (:meth:`TraceFidelity.synthesize`).
   The result is a trace statistically indistinguishable from the measured one
   to the decoder, whose flip times are known exactly.

3. **Score.** Reconstruct each surrogate blind, with the same settings used on
   the measured trace, and score it (:mod:`.analysis`). Repeating over seeds
   gives efficiency, purity, F1 and timing accuracy *with error bars*, at the
   fidelity of the data actually taken.

The efficiency and purity that come out are also what turns the measured flip
count into a rate: a reconstruction reporting ``N`` flips at efficiency ``eps``
and purity ``rho`` implies ``N * rho / eps`` real ones, which is
:attr:`BenchmarkReport.corrected_rate_hz`.

What the surrogate does **not** carry
-------------------------------------
The replay is faithful to what the decoder consumes, and no more. Three limits
are worth stating plainly, because each makes the benchmark *optimistic*:

* **Noise is white and Gaussian**, isotropic in the two quadratures, at the
  sigma read off the trace's minor principal axis. Real amplifier drift, 1/f
  gain wander or interference sit outside that model.
* **The nuisances the pipeline already removed are replayed as removed-able.**
  A ramp-reset comb found on the measured trace is rebuilt into the surrogate,
  so the benchmark does test rediscovering it; a comb the pipeline *missed*
  cannot be, and its cost is invisible here. :attr:`TraceFidelity.reset_comb`
  says which case you are in.
* **The tunnelling process is a pure telegraph** at one rate unless bursts are
  passed explicitly. Burst crowding is the dominant efficiency loss when it is
  present, so pass a :class:`~qpd.simulator.QuasiparticleBurstModel` when the
  data has bursts in it -- see ``bursts=`` on :func:`benchmark_reconstruction`.

A benchmark run on a trace flagged ``degenerate`` is meaningless in a specific
and dangerous way: the fitted model is spurious, so the surrogate replays a
*fiction* which the decoder then reconstructs perfectly. The report carries the
flag through and :meth:`BenchmarkReport.summary` says so.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from qpd.simulator.parity import EVEN, generate_parity_trajectory

from .analysis import DEFAULT_SCALE, DEFAULT_TOL, FlipScore, score_flips
from .bursts import detect_bursts, match_bursts
from .emission import EmissionModel, learn_emission_model, validate_trace
from .reconstruct import reconstruct_parity_flips_ramped
from .static_bias import StaticBlobModel, reconstruct_parity_flips_static

__all__ = [
    "as_complex_trace",
    "TraceFidelity",
    "characterize_trace",
    "BenchmarkTrial",
    "BenchmarkReport",
    "benchmark_reconstruction",
    "benchmark_vs_noise",
    "sweep_rate",
    "sweep_burst_size",
    "BurstSizePoint",
    "implied_rate_hz",
    "burst_n50",
]

# Fractional departure of `BenchmarkReport.closure` from 1 that earns a warning.
# Measured on the reference device: 1.5% at a usable fixed bias (n_g = 0.20),
# rising to 10% once the contrast falls near 1 and the mixture fit starts
# overstating the separation.
_CLOSURE_TOL = 0.05


def as_complex_trace(data) -> np.ndarray:
    """Coerce a user-supplied readout trace to a 1-D complex array.

    Accepts what measured data actually arrives as: an already-complex 1-D
    array, or a real array of I and Q laid out as ``(n, 2)`` or ``(2, n)``.
    A real 1-D array is rejected rather than silently treated as ``Q = 0`` --
    that would make the trace lie on a line, and a phase-only measurement has
    no business going through this pipeline unannounced.
    """
    arr = np.asarray(data)
    if np.iscomplexobj(arr):
        if arr.ndim != 1:
            raise ValueError(f"complex trace must be 1-D; got shape {arr.shape}")
        return arr.astype(complex)
    if arr.ndim == 2:
        if arr.shape[1] == 2:
            return arr[:, 0].astype(float) + 1j * arr[:, 1].astype(float)
        if arr.shape[0] == 2:
            return arr[0].astype(float) + 1j * arr[1].astype(float)
        raise ValueError(
            f"real 2-D input must be (n, 2) or (2, n) I/Q; got {arr.shape}")
    raise ValueError(
        "trace must be complex 1-D, or real (n, 2) / (2, n) I and Q; got a "
        f"real array of shape {arr.shape}. If you have I and Q separately, "
        "pass `i + 1j*q`.")


@dataclass
class TraceFidelity:
    """The measured trace's fidelity, in the form the decoder sees it.

    This is the fitted forward model plus the numbers that set how hard the
    reconstruction is: the branch separation in units of the noise, the
    tunnelling rate, and (when ``n_g`` is swept) the ramp structure. It is what
    :meth:`synthesize` replays.
    """

    mode: str  # "static" (fixed n_g) or "ramped" (swept n_g)
    sample_rate: float
    n_samples: int
    t0: float
    sigma: float  # per-quadrature noise std, in trace units
    contrast_median: float  # branch separation / sigma (median over the trace)
    contrast_max: float
    rate_hz: float  # decoded per-state tunnelling rate
    n_flips: int  # flips found on the measured trace
    degenerate: bool
    sample_fidelity: float
    decoded_fidelity: float
    fold_period: float | None = None  # [s] swept case only
    reset_comb: object | None = None  # ResetComb found on the measured trace
    charge_jump_times: np.ndarray = field(
        default_factory=lambda: np.empty(0, float))
    model: EmissionModel | StaticBlobModel | None = None
    result: object | None = None  # the reconstruction of the measured trace
    recon_kwargs: dict = field(default_factory=dict)

    @property
    def duration(self) -> float:
        return self.n_samples / self.sample_rate

    def _branch_means(self, n: int) -> tuple[np.ndarray, np.ndarray]:
        """Per-sample means of the two branches on the projected axis."""
        if self.mode == "static":
            m = self.model
            return np.full(n, m.mu_a), np.full(n, m.mu_b)
        a, b = self.model.branch_means()
        if a.size < n:
            raise ValueError(
                f"cannot synthesize {n} samples from a model fitted on "
                f"{a.size}; the ramp profile is only defined where it was "
                "learned. Use n_samples <= the measured length.")
        return a[:n], b[:n]

    def synthesize(
        self,
        seed: int | None = None,
        *,
        n_samples: int | None = None,
        rate_hz: float | None = None,
        noise_scale: float = 1.0,
        bursts=None,
        return_bursts: bool = False,
    ):
        """Draw one surrogate trace carrying this fidelity, with known truth.

        A fresh telegraph trajectory is pushed through the *fitted* emission
        model and given fresh noise: the branch means, the splitting profile,
        the ramp reset schedule and the noise sigma are all the measured ones,
        so the decoder faces the same discrimination problem it faced on the
        real data.

        The noise is isotropic in the plane at the fitted ``sigma``, which is
        the model the pipeline itself assumes -- the branch structure lives
        along one axis and the perpendicular axis is noise only. Off-axis
        structure in the measured trace is therefore not reproduced.

        Parameters
        ----------
        seed : int, optional
            Seed for the trajectory and the noise.
        n_samples : int, optional
            Length of the surrogate; defaults to the measured length. May be
            shorter (a prefix), which is the way to keep a benchmark cheap on a
            long trace, but not longer -- the ramp profile is only defined over
            the window it was learned on.
        rate_hz : float, optional
            Per-state tunnelling rate to inject. Defaults to the rate decoded
            from the measured trace.
        noise_scale : float
            Multiply the fitted sigma. This is the "what if the readout were
            better" knob: ``0.5`` halves the noise and doubles every contrast.
        bursts : QuasiparticleBurstModel, optional
            Superimpose quasiparticle bursts on the background telegraph. Pass
            one when the measured data contains bursts -- crowded flips are the
            dominant efficiency loss and the plain telegraph will not show it.
        return_bursts : bool
            Also return the per-burst :class:`~qpd.simulator.BurstTruth`
            records, which carry each burst's true multiplicity and extent --
            the truth a burst-level study is scored against.

        Returns
        -------
        (iq, flip_times) : tuple of arrays
            The surrogate trace, and the exact sub-sample truth flip times on
            the same time axis as the measured trace (offset by ``t0``).
            With ``return_bursts=True``, ``(iq, flip_times, burst_truth)``.
        """
        if self.model is None:
            raise ValueError("no fitted model to replay")
        n = int(self.n_samples if n_samples is None else n_samples)
        if n < 8:
            raise ValueError("need at least 8 samples")
        rng = np.random.default_rng(seed)
        dt = 1.0 / self.sample_rate
        t = np.arange(n) * dt

        rate = float(self.rate_hz if rate_hz is None else rate_hz)
        if not np.isfinite(rate) or rate < 0:
            raise ValueError(f"rate_hz must be finite and non-negative; got {rate}")

        extra, burst_truth = None, []
        if bursts is not None:
            extra, burst_truth = bursts.sample(rng)
        parity, flips = generate_parity_trajectory(
            t, rate, rate, rng, extra_flip_times=extra, return_flip_times=True)

        mu_a, mu_b = self._branch_means(n)
        x = np.where(parity == EVEN, mu_a, mu_b)

        sigma = float(self.sigma) * float(noise_scale)
        noise = rng.normal(0.0, sigma, n) + 1j * rng.normal(0.0, sigma, n)
        iq = self.model.origin + self.model.direction * (x + noise)
        if return_bursts:
            return iq, flips + self.t0, burst_truth
        return iq, flips + self.t0

    def describe(self) -> str:
        """One-block human summary of the measured trace's fidelity."""
        lines = [
            f"mode                 {self.mode}",
            f"samples              {self.n_samples}  "
            f"({self.duration:.4g} s at {self.sample_rate:.4g} Hz)",
            f"noise sigma          {self.sigma:.4g}",
            f"contrast (median)    {self.contrast_median:.3f}",
            f"contrast (max)       {self.contrast_max:.3f}",
            f"decoded rate         {self.rate_hz:.4g} Hz per state",
            f"flips found          {self.n_flips}",
            f"single-sample fid.   {self.sample_fidelity:.4f}",
            f"decoded fidelity     {self.decoded_fidelity:.4f}",
        ]
        if self.mode == "ramped":
            lines.append(f"fold period          {self.fold_period:.6g} s")
            lines.append(
                "ramp reset comb      "
                + (f"period {self.reset_comb.period:.6g} s"
                   if self.reset_comb is not None else "none found"))
            lines.append(f"charge jumps         {self.charge_jump_times.size}")
        if self.degenerate:
            lines.append("DEGENERATE           yes -- the fitted model is not "
                         "trustworthy; see below")
        return "\n".join(lines)


def characterize_trace(
    iq,
    sample_rate: float,
    *,
    mode: str = "auto",
    t0: float = 0.0,
    min_ramp_cycles: float = 50.0,
    **recon_kwargs,
) -> TraceFidelity:
    """Fit the measured trace and package what sets the reconstruction's fidelity.

    Runs the blind reconstruction (which is what learns the forward model) and
    returns it alongside the numbers the surrogate needs. The reconstruction of
    the measured trace is kept on :attr:`TraceFidelity.result`, so this is not a
    wasted pass -- it is the analysis of the real data.

    Parameters
    ----------
    iq : array
        The measured trace. Complex 1-D, or real ``(n, 2)`` / ``(2, n)`` I/Q
        (see :func:`as_complex_trace`).
    sample_rate : float
        Sampling rate [Hz].
    mode : {"auto", "static", "ramped"}
        Which entry point applies. ``"auto"`` decides by how many fold cycles
        the trace holds (see ``min_ramp_cycles``). Override it when you know
        how the measurement was driven -- that is always the safer choice.
    t0 : float
        Time of the first sample [s].
    min_ramp_cycles : float
        Fold cycles the trace must contain before ``mode="auto"`` calls it
        swept. The presence of a spectral peak is *not* sufficient: a fixed
        bias produces one at the search band's lower edge from the telegraph's
        own Lorentzian spectrum, and following it sends the trace down the
        ramped pipeline and inflates the decoded rate several-fold. A real
        sweep packs thousands of cycles into a trace, so this threshold sits
        two orders of magnitude clear of both cases.
    **recon_kwargs
        Passed to the reconstruction, and *reused unchanged on every surrogate*
        so the benchmark measures the settings you will actually run. E.g.
        ``ramp_period=2e-3`` (swept) or ``segment_blocks=8`` (fixed bias).
    """
    iq = validate_trace(as_complex_trace(iq), sample_rate)
    dt = 1.0 / float(sample_rate)

    if mode == "auto":
        # A spectral peak alone does NOT mean the bias was swept. The telegraph
        # itself has a Lorentzian spectrum whose in-band maximum sits at the
        # bottom edge of the period search, and on a fixed-bias trace that is
        # reported as a highly "significant" period of order the trace length
        # (measured: 0.87 s on a 5 s trace, prominence 118). Taking it at face
        # value sends a fixed-bias trace down the ramped pipeline, which then
        # models telegraph noise as a ramp and returns a rate several times too
        # high -- silently.
        #
        # A real n_g sweep is distinguished by *how many* fold cycles it packs
        # into the trace: thousands (27,500 at the reference 500 Hz ramp, 2,750
        # at 50 Hz), against the ~5 the search band's lower edge allows. The
        # cycle count separates the two by two orders of magnitude either side
        # of this threshold.
        probe = learn_emission_model(iq, dt)
        cycles = (iq.size * dt / probe.fold_period
                  if probe.fold_period else 0.0)
        mode = "ramped" if cycles >= float(min_ramp_cycles) else "static"
    if mode not in ("static", "ramped"):
        raise ValueError(f"mode must be 'auto', 'static' or 'ramped'; got {mode!r}")

    if mode == "static":
        res = reconstruct_parity_flips_static(iq, sample_rate, t0=t0,
                                              **recon_kwargs)
        contrast_med = contrast_max = float(res.contrast)
        model = res.model
        fold_period = None
        comb = None
        jumps = np.empty(0, float)
        sigma = float(model.sigma)
    else:
        res = reconstruct_parity_flips_ramped(iq, sample_rate, t0=t0,
                                              **recon_kwargs)
        c = np.asarray(res.contrast, dtype=float)
        contrast_med = float(np.median(c)) if c.size else float("nan")
        contrast_max = float(np.max(c)) if c.size else float("nan")
        model = res.emission
        fold_period = model.fold_period
        comb = res.reset_comb
        jumps = res.charge_jump_times
        sigma = float(model.sigma)

    return TraceFidelity(
        mode=mode,
        sample_rate=float(sample_rate),
        n_samples=int(iq.size),
        t0=float(t0),
        sigma=sigma,
        contrast_median=contrast_med,
        contrast_max=contrast_max,
        rate_hz=float(res.rate_hz),
        n_flips=int(res.flip_times.size),
        degenerate=bool(res.degenerate),
        sample_fidelity=float(res.sample_fidelity),
        decoded_fidelity=float(res.decoded_fidelity),
        fold_period=fold_period,
        reset_comb=comb,
        charge_jump_times=np.asarray(jumps, dtype=float),
        model=model,
        result=res,
        recon_kwargs=dict(recon_kwargs),
    )


@dataclass
class BenchmarkTrial:
    """One surrogate: its score against truth, and its own self-diagnostics."""

    seed: int
    score: FlipScore
    degenerate: bool
    rate_injected: float  # rate this surrogate was drawn at
    rate_hz: float  # rate decoded back out of the surrogate
    contrast_median: float  # contrast re-measured on the surrogate
    n_truth_injected: int


@dataclass
class BenchmarkReport:
    """Reconstruction performance at the measured trace's fidelity.

    The aggregates are over surrogate trials, so the spread is the run-to-run
    scatter at this fidelity -- what an error bar on a single measured trace's
    efficiency should be.
    """

    fidelity: TraceFidelity
    trials: list[BenchmarkTrial]
    tol: float = DEFAULT_TOL
    scale: float = DEFAULT_SCALE
    noise_scale: float = 1.0
    rate_jitter: bool = True
    select_min_contrast: float | None = None
    warnings: list[str] = field(default_factory=list)

    @property
    def selected(self) -> list:
        """Trials that pass ``select_min_contrast``, i.e. the ones scored.

        A surrogate is subject to the same acceptance cut you would apply to a
        real chunk. That matters most at low injected rate, where a short trace
        can contain *no* flips at all: with no transition the two-branch model
        is unidentifiable, EM splits a single blob spuriously (contrast ~0.9),
        and the decoder segments noise into hundreds of fabricated events.
        Pooling those into the average is what makes purity appear to collapse
        at low rate -- measured on a 1 s, 10 kSa/s trace at contrast 2.4, 36%
        of 1 Hz surrogates held no flips and purity read 0.08; applying the
        same contrast > 1.7 cut used to select the real chunks gives 1.00.
        """
        if self.select_min_contrast is None:
            return list(self.trials)
        thr = float(self.select_min_contrast)
        return [t for t in self.trials
                if np.isfinite(t.contrast_median) and t.contrast_median > thr]

    @property
    def n_excluded(self) -> int:
        """Trials rejected by the acceptance cut."""
        return len(self.trials) - len(self.selected)

    @property
    def selection_fraction(self) -> float:
        """Fraction of surrogates that pass the cut.

        This is a *result*, not bookkeeping: at a given rate and trace length it
        is the probability that a real chunk will be usable at all.
        """
        return len(self.selected) / len(self.trials) if self.trials else np.nan

    def _agg(self, attr: str) -> tuple[float, float]:
        v = np.array([getattr(t.score, attr) for t in self.selected], dtype=float)
        v = v[np.isfinite(v)]
        if v.size == 0:
            return float("nan"), float("nan")
        return float(np.mean(v)), float(np.std(v, ddof=1) if v.size > 1 else 0.0)

    def _pooled(self, denom: str) -> tuple[float, float]:
        """Ratio of summed counts, with its binomial error.

        Every trial contributes its *events* rather than its ratio, which
        changes two things and both matter:

        * The error shrinks as ``1/sqrt(total events)``, so it falls with
          ``n_trials`` -- unlike :meth:`_agg`, whose spread converges to the
          population scatter and stays put however long you run.
        * Trials are weighted by the evidence they carry. Averaging ratios
          gives a trial holding 3 events the same weight as one holding 500,
          which hides rare catastrophic trials. Measured at 1 Hz on a
          contrast-1.6 trace: two surrogates out of 200 contained *no* real
          flips at all, yet the decoder segmented noise into ~100 spurious ones
          each. They are 1% of the mean-of-ratios and 26% of every prediction
          made -- mean-of-ratios reported purity 0.96, pooled reported 0.74.

        The interval is the normal approximation, so it is only sensible away
        from 0 and 1; near the bounds prefer Wilson or Clopper-Pearson.
        """
        num = sum(int(t.score.n_matched) for t in self.selected)
        den = sum(int(getattr(t.score, denom)) for t in self.selected)
        if den <= 0:
            return float("nan"), float("nan")
        p = num / den
        return float(p), float(np.sqrt(max(p * (1.0 - p), 0.0) / den))

    def _ratio_clustered(self, denom: str) -> tuple[float, float]:
        """Pooled ratio with the linearised cluster variance -- closed form.

        The delta-method variance of a ratio of cluster totals (Cochran's
        ratio estimator, standard in survey sampling):

            V(p) = n / ((n-1) * (sum b_i)^2) * sum_i (a_i - p*b_i)^2

        with ``a_i`` the matched count in trial ``i`` and ``b_i`` its
        denominator. Each trial enters as one *cluster*, so this makes the same
        correction as :meth:`_bootstrap` -- the effective sample size is the
        number of trials, not of events -- but it is deterministic, exact and
        instant: no resample count to pick, no seed.

        The two agree to 3-4 significant figures wherever the trial sizes are
        not wildly heterogeneous (measured: 0.00305 vs 0.00304 for purity at
        3 Hz, 0.00124 vs 0.00123 at 100 Hz). They part company exactly where
        the linearisation is least trustworthy -- at 1 Hz, where two trials of
        200 carry a quarter of all events, delta gives 0.114, the jackknife
        0.126 and the bootstrap 0.112. Treat that spread as the warning it is,
        and cross-check with :attr:`purity_bootstrap` when it appears.
        """
        a = np.array([t.score.n_matched for t in self.selected], dtype=float)
        b = np.array([getattr(t.score, denom) for t in self.selected],
                     dtype=float)
        n = a.size
        tot = b.sum()
        if n == 0 or tot <= 0:
            return float("nan"), float("nan")
        p = float(a.sum() / tot)
        if n < 2:
            return p, float("nan")
        var = n / ((n - 1) * tot ** 2) * float(np.sum((a - p * b) ** 2))
        return p, float(np.sqrt(max(var, 0.0)))

    @property
    def efficiency_clustered(self) -> tuple[float, float]:
        """Pooled recall with the closed-form cluster error. **The default.**

        Estimates the same thing as :attr:`efficiency_bootstrap` without a
        resample count or a seed.
        """
        return self._ratio_clustered("n_truth")

    @property
    def purity_clustered(self) -> tuple[float, float]:
        """Pooled precision with the closed-form cluster error."""
        return self._ratio_clustered("n_pred")

    # -- detection ---------------------------------------------------------
    @property
    def efficiency(self) -> tuple[float, float]:
        """Recall (matched / true) as ``(mean, std)`` over trials.

        The second entry is the trial-to-trial *scatter*, not the precision of
        the mean: it answers "how much would one more trace like mine vary",
        and it does **not** shrink with ``n_trials``. For "how efficient is the
        algorithm here", which does improve with more trials, use
        :attr:`efficiency_clustered`.
        """
        return self._agg("efficiency")

    @property
    def purity(self) -> tuple[float, float]:
        """Precision (matched / predicted) as ``(mean, std)``.

        Same caveat as :attr:`efficiency`; see :attr:`purity_pooled`.
        """
        return self._agg("purity")

    def _bootstrap(self, denom: str, n_resamples: int = 2000,
                   seed: int = 0) -> tuple[float, float]:
        """Pooled ratio with a cluster-aware error, by resampling *trials*.

        The binomial error of :meth:`_pooled` assumes the pooled events are
        independent. They are not: they arrive in trial-sized clusters, and a
        single surrogate that segments noise contributes a hundred correlated
        failures at once. The effective sample size is therefore the number of
        *trials*, not the number of events, and the binomial interval is far
        too tight -- measured at 1 Hz, the pooled purity wandered over
        0.98 / 0.77 / 0.83 for 25 / 100 / 400 trials while the binomial error
        claimed +/-0.01.

        Resampling whole trials with replacement respects that clustering, and
        still tightens as ``1/sqrt(n_trials)``. It costs nothing: the per-trial
        counts are already stored, so no reconstruction is repeated.
        """
        num = np.array([t.score.n_matched for t in self.selected], dtype=float)
        den = np.array([getattr(t.score, denom) for t in self.selected],
                       dtype=float)
        tot = den.sum()
        if tot <= 0 or num.size == 0:
            return float("nan"), float("nan")
        point = float(num.sum() / tot)
        rng = np.random.default_rng(seed)
        idx = rng.integers(0, num.size, size=(int(n_resamples), num.size))
        d = den[idx].sum(axis=1)
        ok = d > 0
        if not np.any(ok):
            return point, float("nan")
        draws = num[idx].sum(axis=1)[ok] / d[ok]
        return point, float(np.std(draws, ddof=1))

    @property
    def efficiency_bootstrap(self) -> tuple[float, float]:
        """Pooled recall with a trial-level bootstrap error.

        The most defensible of the three: it improves with ``n_trials`` like
        :attr:`efficiency_pooled`, but unlike it does not pretend the events
        within a trial are independent.
        """
        return self._bootstrap("n_truth")

    @property
    def purity_bootstrap(self) -> tuple[float, float]:
        """Pooled precision with a trial-level bootstrap error."""
        return self._bootstrap("n_pred")

    @property
    def efficiency_pooled(self) -> tuple[float, float]:
        """Recall over pooled counts, ``(value, binomial error)``.

        Use this when the question is *how efficient is the reconstruction at
        this operating point* -- a property of the algorithm, estimated by
        Monte Carlo, where more trials should buy a better answer. Use
        :attr:`efficiency` when the question is *how much would one more trace
        like mine vary*.
        """
        return self._pooled("n_truth")

    @property
    def purity_pooled(self) -> tuple[float, float]:
        """Precision over pooled counts, ``(value, binomial error)``."""
        return self._pooled("n_pred")

    @property
    def hard_f1_pooled(self) -> float:
        """F1 built from the pooled precision and recall."""
        e, _ = self.efficiency_pooled
        p, _ = self.purity_pooled
        if not (np.isfinite(e) and np.isfinite(p)) or (e + p) <= 0:
            return float("nan")
        return float(2 * e * p / (e + p))

    @property
    def hard_f1(self) -> tuple[float, float]:
        return self._agg("hard_f1")

    @property
    def soft_f1(self) -> tuple[float, float]:
        """F1 with each match weighted by its timing residual."""
        return self._agg("soft_f1")

    # -- timing ------------------------------------------------------------
    @property
    def timing_bias_s(self) -> tuple[float, float]:
        """Mean signed residual ``pred - truth``, as ``(mean, std)`` [s]."""
        return self._agg("bias_s")

    @property
    def timing_rms_s(self) -> tuple[float, float]:
        return self._agg("rms_s")

    @property
    def residuals(self) -> np.ndarray:
        """All matched timing residuals, pooled over trials [s]."""
        parts = [t.score.dt for t in self.selected if t.score.dt.size]
        return np.concatenate(parts) if parts else np.empty(0, float)

    # -- what this implies for the measured trace --------------------------
    @property
    def corrected_rate_hz(self) -> float:
        """Tunnelling rate implied by the measured flip count, de-biased.

        The decoder misses flips (efficiency < 1) and invents them (purity < 1),
        and the two do not cancel. Since
        ``n_pred * purity / efficiency == n_truth`` identically, applying the
        benchmark's efficiency and purity to the flip count the pipeline
        actually reported inverts that bias.

        The count has to be the one the pipeline *output*, not
        :attr:`TraceFidelity.rate_hz`. Those differ whenever a post-decode
        filter is in play: ``rate_hz`` comes from the HMM's ``p_flip``, which
        is estimated *before* ``min_confidence`` drops any flips, while
        efficiency and purity are measured on surrogates with the filter
        applied. Mixing them overstates the rate by the fraction the cut
        removes -- measured on the reference device at ``min_confidence=0.4``,
        that was +35% to +52% across six traces (mean absolute error 42%,
        against 4% using the reported count).

        The pooled efficiency and purity are used rather than the
        mean-of-ratios, since this is an inversion of ratios of totals.
        """
        eff, _ = self.efficiency_clustered
        pur, _ = self.purity_clustered
        duration = self.fidelity.duration
        if (not (np.isfinite(eff) and np.isfinite(pur))
                or eff <= 0 or duration <= 0):
            return float("nan")
        return float(self.fidelity.n_flips / duration * pur / eff)

    @property
    def closure(self) -> float:
        """Surrogate contrast / measured contrast; should be ~1.

        The one check that the replay actually reproduced the measurement. A
        value far from 1 means the surrogate is not at the data's fidelity and
        the numbers above do not describe it.

        **Defined only for** ``noise_scale == 1``, and ``nan`` otherwise. The
        test works because both sides are measured by the same estimator at the
        same fidelity, so its bias cancels. Under a deliberate noise rescaling
        they sit at *different* contrasts, and the mixture fit is not linear
        between them -- below a true contrast of about 1 it splits a single
        blob into a spurious pair and floors out near 0.93 (measured: a
        surrogate at a true 0.59 reads back 0.93). Comparing across that would
        report a replay failure where the replay is exact; the injected noise
        itself scales to within 0.3%.
        """
        if self.noise_scale != 1.0:
            return float("nan")
        c = np.array([t.contrast_median for t in self.selected], dtype=float)
        c = c[np.isfinite(c)]
        ref = self.fidelity.contrast_median
        if c.size == 0 or not np.isfinite(ref) or ref == 0:
            return float("nan")
        return float(np.mean(c) / ref)

    def summary(self) -> str:
        """Printable report: the measured fidelity, then performance at it."""
        f = self.fidelity
        eff, eff_e = self.efficiency
        pur, pur_e = self.purity
        h, h_e = self.hard_f1
        s, s_e = self.soft_f1
        b, b_e = self.timing_bias_s
        r, r_e = self.timing_rms_s
        n_deg = sum(t.degenerate for t in self.trials)
        out = [
            "Measured trace",
            "--------------",
            f.describe(),
            "",
            f"Reconstruction performance at this fidelity "
            f"({len(self.trials)} surrogate trials"
            + (f", noise x{self.noise_scale:g}" if self.noise_scale != 1.0 else "")
            + ("" if self.rate_jitter else ", rate pinned")
            + ")",
            "-" * 72,
            f"efficiency (recall)  {eff:.3f} +/- {eff_e:.3f}",
            f"purity (precision)   {pur:.3f} +/- {pur_e:.3f}",
            f"hard F1              {h:.3f} +/- {h_e:.3f}",
            f"soft F1 (timed)      {s:.3f} +/- {s_e:.3f}",
            f"timing bias          {b * 1e6:+.1f} +/- {b_e * 1e6:.1f} us",
            f"timing rms           {r * 1e6:.1f} +/- {r_e * 1e6:.1f} us",
            f"match tolerance      {self.tol * 1e6:.0f} us",
            (f"accepted             {len(self.selected)}/{len(self.trials)} "
             f"surrogates (contrast > {self.select_min_contrast:g})"
             if self.select_min_contrast is not None else
             f"accepted             all {len(self.trials)} surrogates"),
            f"closure (sim/meas)   {self.closure:.3f}",
            "",
            f"Implied rate on the measured trace: "
            f"{self.corrected_rate_hz:.4g} Hz per state "
            f"(decoder reported {f.rate_hz:.4g} Hz)",
            "",
            "Per-trial",
            "  " + FlipScore.header(),
        ]
        for t in self.trials:
            out.append(f"  {t.score.row()}")
        if n_deg:
            out.append("")
            out.append(f"NOTE: {n_deg}/{len(self.trials)} surrogates were "
                       "flagged degenerate by the reconstruction.")
        for w in self.warnings:
            out.append(f"WARNING: {w}")
        return "\n".join(out)


def _reconstruct(iq, fidelity: TraceFidelity):
    """Run the same entry point and settings used on the measured trace."""
    if fidelity.mode == "static":
        return reconstruct_parity_flips_static(
            iq, fidelity.sample_rate, t0=fidelity.t0, **fidelity.recon_kwargs)
    return reconstruct_parity_flips_ramped(
        iq, fidelity.sample_rate, t0=fidelity.t0, **fidelity.recon_kwargs)


def benchmark_reconstruction(
    iq=None,
    sample_rate: float | None = None,
    *,
    fidelity: TraceFidelity | None = None,
    n_trials: int = 8,
    seed: int = 0,
    trial_samples: int | None = None,
    trial_duration: float | None = None,
    noise_scale: float = 1.0,
    rate_hz: float | None = None,
    rate_jitter: bool = True,
    select_min_contrast: float | None = None,
    bursts=None,
    tol: float = DEFAULT_TOL,
    scale: float = DEFAULT_SCALE,
    calibrate_rate: bool = False,
    mode: str = "auto",
    t0: float = 0.0,
    **recon_kwargs,
) -> BenchmarkReport:
    """Measure reconstruction efficiency and timing accuracy on *your* data.

    Give it a measured I/Q trace. It fits that trace, replays its fidelity into
    surrogate traces that do have truth, reconstructs those blind with the same
    settings, and scores them -- so the efficiency and accuracy that come back
    describe the reconstruction of the data you actually took, not a generic
    device.

    Parameters
    ----------
    iq : array
        The measured trace: complex 1-D, or real ``(n, 2)`` / ``(2, n)`` I/Q.
        Omit only when passing a pre-computed ``fidelity``.
    sample_rate : float
        Sampling rate [Hz].
    fidelity : TraceFidelity, optional
        A characterisation from :func:`characterize_trace`. Pass one to reuse a
        fit across several benchmark settings (e.g. a noise sweep) instead of
        re-reconstructing the measured trace every time.
    n_trials : int
        Number of surrogate traces. The spread across them is the reported
        error bar, so a handful is enough for a mean and about a dozen before
        the error bar itself is stable.
    seed : int
        Base seed; trial ``k`` uses ``seed + k``.
    trial_samples, trial_duration : int / float, optional
        Shorten each surrogate to this many samples (or seconds). The swept
        pipeline costs a few seconds per second of trace, so this is the knob
        that keeps a benchmark of a long measurement affordable. Surrogates
        cannot be made *longer* than the measured trace -- the fitted ramp
        profile is only defined over the window it was learned on. Note a short
        surrogate holds fewer flips, which widens the error bars.
    noise_scale : float
        Scale the fitted noise sigma. ``1.0`` benchmarks the measurement as
        taken; other values answer "what would a better (or worse) readout
        buy?" -- see :func:`benchmark_vs_noise`.
    rate_hz : float, optional
        Override the injected tunnelling rate. Defaults to the rate decoded
        from the measured trace.
    rate_jitter : bool
        Draw each surrogate's rate from the Poisson posterior for the flip
        count actually seen, ``Gamma(N + 1/2, 1/T)``, instead of pinning every
        one at the decoded rate. This is on by default because the decoded rate
        is an ``N``-event estimate: a trace holding a dozen flips knows its own
        rate only to ~30%, and how often two flips crowd into one dwell -- the
        dominant efficiency loss -- is very sensitive to it. Pinning the rate
        produces error bars that describe one particular realisation rather
        than the measurement, and reports things like ``1.000 +/- 0.000`` off a
        nine-flip trace. Turn it off to isolate the decoder's own scatter at a
        fixed rate.
    select_min_contrast : float, optional
        Accept a surrogate only if its *re-fitted* contrast exceeds this,
        applying to the simulated traces the same cut you would apply when
        choosing which measured chunks to analyse. Excluded trials are counted
        (:attr:`BenchmarkReport.n_excluded`) and reported, never silently
        dropped, because the rejection rate is itself a result -- at a given
        rate and trace length it is the probability that a real chunk is usable
        at all.

        Strongly recommended whenever the sweep reaches low rates. A trace short
        enough to contain no flips has no identifiable two-branch model, and the
        decoder fabricates hundreds of events from noise; pooling those makes
        purity appear to collapse at low rate when the reconstruction is fine.
    bursts : QuasiparticleBurstModel, optional
        Superimpose quasiparticle bursts on each surrogate. Pass one if the
        measurement contains bursts: crowded flips dominate the efficiency loss
        and a plain telegraph will not reveal it.
    tol : float
        Hard matching tolerance [s]; a prediction further than this from a true
        flip is a miss *and* a false positive. The default is the grader's.
    scale : float
        Timing kernel scale [s] for the soft F1.
    calibrate_rate : bool
        Run one extra round of trials first and use its efficiency and purity
        to de-bias the injected rate. Worth it when the decoded rate is
        suspect, which is exactly when the contrast is marginal; it doubles the
        cost. It is skipped outright on a degenerate trace, where the decoded
        "rate" is not a flip count and the correction diverges.
    mode, t0, **recon_kwargs
        Passed to :func:`characterize_trace` (and, for ``recon_kwargs``, to
        every surrogate reconstruction as well).

    Returns
    -------
    BenchmarkReport
    """
    if fidelity is None:
        if iq is None or sample_rate is None:
            raise ValueError("pass a measured trace (iq, sample_rate), or a "
                             "pre-computed fidelity=")
        fidelity = characterize_trace(iq, sample_rate, mode=mode, t0=t0,
                                      **recon_kwargs)
    elif recon_kwargs:
        # A pre-computed fidelity already carries its own settings, and the
        # surrogates must be reconstructed with those. Silently dropping a
        # second set here would report the wrong pipeline's performance.
        raise ValueError(
            f"reconstruction settings {sorted(recon_kwargs)} cannot be given "
            "alongside fidelity=; they belong to characterize_trace, and the "
            "supplied fidelity already carries "
            f"{sorted(fidelity.recon_kwargs) or 'none'}")

    n = fidelity.n_samples
    if trial_samples is not None and trial_duration is not None:
        raise ValueError("pass trial_samples or trial_duration, not both")
    if trial_duration is not None:
        trial_samples = int(round(trial_duration * fidelity.sample_rate))
    if trial_samples is not None:
        if trial_samples > n:
            raise ValueError(
                f"trial_samples={trial_samples} exceeds the measured length "
                f"{n}; the fitted ramp profile is only defined there")
        n = int(trial_samples)

    warnings: list[str] = []
    if fidelity.degenerate:
        warnings.append(
            "the measured trace is flagged degenerate: the fitted model is "
            "spurious, so these surrogates replay a fiction and the scores "
            "below describe that fiction, not the measurement.")
    if fidelity.mode == "ramped" and fidelity.reset_comb is None:
        warnings.append(
            "no ramp-reset comb was found on the measured trace. If the sweep "
            "does span an odd number of half Cooper pairs, the resets are "
            "missing from the surrogate too and the benchmark is optimistic; "
            "supply ramp_period= if you know it.")
    if fidelity.mode == "ramped" and fidelity.charge_jump_times.size:
        warnings.append(
            f"{fidelity.charge_jump_times.size} offset-charge jump(s) were "
            "found on the measured trace. Their phase realignment is baked "
            "into the replayed model, so the surrogates do not test "
            "rediscovering them.")

    def _rates(rate, base_seed):
        """One injected rate per trial, jittered by the counting uncertainty."""
        if not rate_jitter or fidelity.n_flips <= 0 or rate <= 0:
            return np.full(int(n_trials), float(rate))
        # Jeffreys posterior for a Poisson rate given N counts in T seconds,
        # rescaled so its mean is the rate actually being injected (which
        # `calibrate_rate` may have moved off the raw decoded value).
        rng = np.random.default_rng(int(base_seed) + 7717)
        shape = fidelity.n_flips + 0.5
        return rate * rng.gamma(shape, 1.0 / shape, int(n_trials))

    def _run(rate, base_seed):
        out = []
        for k, r_k in enumerate(_rates(rate, base_seed)):
            s = int(base_seed) + k
            sim_iq, truth = fidelity.synthesize(
                seed=s, n_samples=n, rate_hz=float(r_k),
                noise_scale=noise_scale, bursts=bursts)
            res = _reconstruct(sim_iq, fidelity)
            c = np.asarray(res.contrast, dtype=float)
            c_med = float(np.median(c)) if c.size else float(res.contrast)
            out.append(BenchmarkTrial(
                seed=s,
                score=score_flips(truth, res.flip_times, tol=tol, scale=scale),
                degenerate=bool(res.degenerate),
                rate_injected=float(r_k),
                rate_hz=float(res.rate_hz),
                contrast_median=c_med,
                n_truth_injected=int(truth.size),
            ))
        return out

    rate = fidelity.rate_hz if rate_hz is None else float(rate_hz)
    if calibrate_rate and fidelity.degenerate:
        # The correction is `decoded * purity / efficiency`, which assumes the
        # decoded rate counts flips. On a degenerate trace it counts noise
        # segmentations, the surrogate replays *those* as truth, and the
        # iteration runs away from the device rate rather than toward it
        # (measured at the parity-blind bias: 3.8 kHz "corrected" to 8.0 kHz
        # against a true 11 Hz). Refuse rather than report that number.
        warnings.append(
            "rate calibration skipped: the measured trace is degenerate, so "
            "the decoded rate is not a flip count and correcting it diverges.")
    elif calibrate_rate:
        pilot = _run(rate, seed + 10_000)
        e = np.array([t.score.efficiency for t in pilot], float)
        p = np.array([t.score.purity for t in pilot], float)
        e, p = e[np.isfinite(e)], p[np.isfinite(p)]
        if e.size and p.size and e.mean() > 0:
            corrected = rate * float(p.mean()) / float(e.mean())
            warnings.append(
                f"rate calibrated from {rate:.4g} Hz to {corrected:.4g} Hz "
                "using a pilot round.")
            rate = corrected

    report = BenchmarkReport(
        fidelity=fidelity, trials=_run(rate, seed), tol=tol, scale=scale,
        noise_scale=noise_scale, rate_jitter=bool(rate_jitter),
        select_min_contrast=select_min_contrast, warnings=warnings,
    )
    if report.n_excluded:
        report.warnings.append(
            f"{report.n_excluded}/{len(report.trials)} surrogates were "
            f"rejected by select_min_contrast={select_min_contrast:g} and are "
            f"not in the scores below. That rejection rate is a result in its "
            f"own right: {100 * (1 - report.selection_fraction):.0f}% of "
            f"traces this length at this rate would be unusable.")
    if not report.selected:
        report.warnings.append(
            "no surrogate passed select_min_contrast, so every score is nan.")
    # Closure is the one test that the replay reproduced the measurement, and
    # it is only computable once the trials are in. It drifts above 1 at
    # marginal contrast, where the mixture fit sits on the optimistic side of
    # the true separation and the surrogate then inherits that bias -- so the
    # scores would describe a slightly easier measurement than the real one.
    c = report.closure
    if np.isfinite(c) and abs(c - 1.0) > _CLOSURE_TOL:
        report.warnings.append(
            f"closure is {c:.3f}, not 1: the surrogates re-measure a contrast "
            f"{abs(c - 1.0) * 100:.0f}% "
            f"{'above' if c > 1 else 'below'} the measured trace's, so they "
            f"are {'easier' if c > 1 else 'harder'} than the data. Expect the "
            "scores below to be correspondingly "
            f"{'optimistic' if c > 1 else 'pessimistic'}.")
    return report


def benchmark_vs_noise(
    iq=None,
    sample_rate: float | None = None,
    *,
    noise_scales=(0.5, 1.0, 2.0),
    fidelity: TraceFidelity | None = None,
    mode: str = "auto",
    t0: float = 0.0,
    **kwargs,
) -> list[BenchmarkReport]:
    """Benchmark the same measurement at several readout noise levels.

    ``noise_scale = 1`` is the measurement as taken; below 1 is a quieter
    readout, above 1 a noisier one. Because every contrast scales as
    ``1 / noise_scale``, this maps out how far the present operating point sits
    from the reconstruction's breakdown -- which is the question worth asking
    before spending effort on the amplifier chain.

    The measured trace is characterised once and reused for every point.
    """
    recon = {k: kwargs.pop(k) for k in list(kwargs) if k in _RECON_KEYS}
    if fidelity is None:
        if iq is None or sample_rate is None:
            raise ValueError("pass a measured trace (iq, sample_rate), or a "
                             "pre-computed fidelity=")
        fidelity = characterize_trace(iq, sample_rate, mode=mode, t0=t0,
                                      **recon)
    elif recon:
        raise ValueError(
            f"reconstruction settings {sorted(recon)} cannot be given "
            "alongside fidelity=; the supplied fidelity already carries "
            f"{sorted(fidelity.recon_kwargs) or 'none'}")
    return [benchmark_reconstruction(fidelity=fidelity, noise_scale=float(s),
                                     **kwargs)
            for s in noise_scales]


def sweep_rate(
    fidelity: TraceFidelity,
    rates_hz,
    *,
    n_trials: int = 6,
    seed: int = 0,
    adaptive_tol: bool = True,
    tol_dwell_fraction: float = 0.25,
    **kwargs,
) -> list[BenchmarkReport]:
    """Benchmark the same measurement across background tunnelling rates.

    The background is a Poisson process, so the rate alone decides how often
    two tunnels land close enough to be unresolvable -- and a pair separated by
    less than the decoder can resolve is not merely mistimed, it is *invisible*
    (two toggles return the parity to where it started). Efficiency therefore
    falls with rate no matter how good the readout is, and this maps out where.

    The rate is pinned rather than jittered here: it is the independent
    variable, so the counting uncertainty of the measured trace has no place in
    it. Everything else -- contrast, noise, ramp structure -- stays at the
    measured trace's fidelity.

    Parameters
    ----------
    adaptive_tol : bool
        Shrink the matching tolerance in step with the rate, to
        ``tol_dwell_fraction / rate``, never above the grader's default. This
        is on by default for two reasons that arrive together. Physically, a
        fixed 0.5 ms tolerance is nonsense once flips arrive every 33 us -- it
        would credit a prediction that is fifteen dwells away as a match.
        Computationally, a tolerance wider than the mean gap fuses the whole
        trace into one connected component of the matching graph, and the
        grader's assignment step goes quadratic on ~10^5 events and effectively
        hangs. Set False only when comparing rates at a genuinely fixed
        tolerance, and then keep the rates low.
    tol_dwell_fraction : float
        Fraction of the mean dwell used as the tolerance when ``adaptive_tol``.

    Returns
    -------
    list of BenchmarkReport, one per rate. Each carries the tolerance it was
    scored at in :attr:`BenchmarkReport.tol`.
    """
    kwargs.setdefault("rate_jitter", False)
    base_tol = kwargs.pop("tol", DEFAULT_TOL)
    # Popped once, outside the loop: popping inside would apply a caller's
    # scale to the first rate only and silently fall back for the rest.
    user_scale = kwargs.pop("scale", None)
    out = []
    for r in np.asarray(rates_hz, dtype=float):
        tol = base_tol
        if adaptive_tol and r > 0:
            tol = min(base_tol, float(tol_dwell_fraction) / float(r))
        out.append(benchmark_reconstruction(
            fidelity=fidelity, rate_hz=float(r), n_trials=n_trials, seed=seed,
            tol=tol, scale=(tol / 2.0 if user_scale is None else user_scale),
            **kwargs))
    return out


@dataclass
class BurstSizePoint:
    """Burst-level performance at one expected quasiparticle multiplicity."""

    n_qp_expected: float  # Poisson mean the bursts were drawn at
    n_bursts: int  # truth bursts that actually injected >= 1 quasiparticle
    n_detected: int
    n_qp_true: np.ndarray  # per matched burst
    n_qp_detected: np.ndarray  # per matched burst
    background_rate_hz: float = float("nan")

    @property
    def efficiency(self) -> float:
        """Fraction of true bursts found -- burst-level, not flip-level."""
        return self.n_detected / self.n_bursts if self.n_bursts else np.nan

    @property
    def efficiency_err(self) -> float:
        """Binomial standard error on :attr:`efficiency`."""
        n = self.n_bursts
        if not n:
            return float("nan")
        p = self.n_detected / n
        return float(np.sqrt(max(p * (1 - p), 0.0) / n))

    @property
    def mean_n_qp_true(self) -> float:
        return float(np.mean(self.n_qp_true)) if self.n_qp_true.size else np.nan

    @property
    def mean_n_qp_detected(self) -> float:
        return (float(np.mean(self.n_qp_detected))
                if self.n_qp_detected.size else np.nan)

    @property
    def bias(self) -> float:
        """Mean detected minus true multiplicity, over *detected* bursts."""
        if not self.n_qp_true.size:
            return float("nan")
        return float(np.mean(self.n_qp_detected - self.n_qp_true))

    @property
    def bias_err(self) -> float:
        if self.n_qp_true.size < 2:
            return float("nan")
        d = self.n_qp_detected - self.n_qp_true
        return float(np.std(d, ddof=1) / np.sqrt(d.size))


def sweep_burst_size(
    fidelity: TraceFidelity,
    n_qp_values,
    *,
    background_rate_hz: float | None = None,
    n_trials: int = 6,
    seed: int = 0,
    burst_spacing: float = 0.25,
    burst_tau: float = 1.0e-3,
    burst_mu: float = 1.2e-3,
    burst_sigma: float = 0.4e-3,
    trial_samples: int | None = None,
    select_min_contrast: float | None = None,
    detect_kwargs: dict | None = None,
) -> list[BurstSizePoint]:
    """Burst-level detection efficiency and multiplicity bias vs burst size.

    This asks a different question from :func:`benchmark_reconstruction`, which
    scores individual flips. Here the object is the *burst*: was the energy
    deposition found at all, and once found, how many of its quasiparticles
    were counted? For a detector those are the numbers that matter.

    At each requested multiplicity, bursts are injected at regular intervals on
    top of the measured trace's background, the surrogate is reconstructed
    blind, :func:`~qpd.reconstruction.bursts.detect_bursts` clusters the
    recovered flips, and the clusters are matched to truth.

    Expect two distinct effects, and read them from the two outputs:

    * **Efficiency turns on** with multiplicity -- a 2-3 quasiparticle burst is
      not statistically distinguishable from a background coincidence, while a
      large one is unmistakable.
    * **Multiplicity saturates.** A big burst crowds its tunnels into a few
      milliseconds, unresolved pairs cancel in the parity, and the recovered
      count falls progressively short of the truth. The bias is what
      :attr:`BurstSizePoint.bias` and the ``n_qp_true`` / ``n_qp_detected``
      pairs record, and it is a property of parity readout rather than of the
      clustering.

    Parameters
    ----------
    fidelity : TraceFidelity
        From :func:`characterize_trace` -- the measured trace's fidelity.
    n_qp_values : sequence of float
        Expected quasiparticle counts (Poisson means) to sweep.
    background_rate_hz : float, optional
        Background tunnelling rate to inject, and the rate the burst finder is
        told to test against. Defaults to the rate decoded from the measured
        trace; pass :attr:`BenchmarkReport.corrected_rate_hz` to use the
        de-biased value instead. Getting this from the data matters -- too low
        a value manufactures bursts from background pairs, too high dissolves
        the small ones.
    n_trials : int
        Surrogate traces per multiplicity. Bursts accumulate across trials, so
        the burst count backing each efficiency is ``n_trials`` times the
        number that fit in one trace.
    burst_spacing : float
        Interval between injected bursts [s]. Must stay well above the burst
        duration (~12 ms on the reference device) so they do not merge.
    burst_tau, burst_mu, burst_sigma : float
        EMG shape of the burst arrival profile. ``burst_tau`` defaults to 1 ms,
        giving a burst that lands inside a few milliseconds; the reference
        device's measured value is 3.7 ms, which spreads one over ~12 ms. The
        two interact with ``detect_bursts``'s linking distance, so change them
        together.
    select_min_contrast : float, optional
        Skip a surrogate whose re-fitted contrast falls below this, exactly as
        in :func:`benchmark_reconstruction`. Its bursts contribute to neither
        the numerator nor the denominator, so the efficiency is conditioned on
        a usable trace rather than diluted by unusable ones.
    detect_kwargs : dict, optional
        Passed to :func:`~qpd.reconstruction.bursts.detect_bursts`, e.g.
        ``{"min_flips": 4}`` to trade efficiency for a lower false-burst rate.

    Returns
    -------
    list of BurstSizePoint, one per entry of ``n_qp_values``.
    """
    from qpd.simulator import QuasiparticleBurstModel

    rate = float(fidelity.rate_hz if background_rate_hz is None
                 else background_rate_hz)
    n = int(fidelity.n_samples if trial_samples is None else trial_samples)
    if n > fidelity.n_samples:
        raise ValueError(
            f"trial_samples={n} exceeds the measured length "
            f"{fidelity.n_samples}")
    duration = n / fidelity.sample_rate
    if burst_spacing <= 0:
        raise ValueError("burst_spacing must be positive")
    # Keep a spacing clear of the trace edges so no burst is truncated.
    onsets = np.arange(burst_spacing, duration - burst_spacing, burst_spacing)
    if onsets.size == 0:
        raise ValueError(
            f"no burst fits in a {duration:.4g} s trace at "
            f"burst_spacing={burst_spacing:.4g} s")
    dk = dict(detect_kwargs or {})

    out: list[BurstSizePoint] = []
    for lam in np.asarray(n_qp_values, dtype=float):
        model = QuasiparticleBurstModel(
            times=onsets, tau=burst_tau, mu=burst_mu, sigma=burst_sigma,
            expected_n_qp=float(lam))
        n_true = n_det = 0
        tr, de = [], []
        for k in range(int(n_trials)):
            sim_iq, _, truth = fidelity.synthesize(
                seed=int(seed) + k, n_samples=n, rate_hz=rate,
                bursts=model, return_bursts=True)
            rec = _reconstruct(sim_iq, fidelity)
            if select_min_contrast is not None:
                c = np.asarray(rec.contrast, dtype=float)
                c_med = float(np.median(c)) if c.size else float(rec.contrast)
                if not (np.isfinite(c_med) and c_med > float(select_min_contrast)):
                    continue
            found = detect_bursts(rec.flip_times - fidelity.t0, rate,
                                  duration=duration, **dk)
            for m in match_bursts(truth, found):
                n_true += 1
                if m.detected:
                    n_det += 1
                    tr.append(m.n_qp_true)
                    de.append(m.n_qp_detected)
        out.append(BurstSizePoint(
            n_qp_expected=float(lam), n_bursts=n_true, n_detected=n_det,
            n_qp_true=np.asarray(tr, dtype=float),
            n_qp_detected=np.asarray(de, dtype=float),
            background_rate_hz=rate,
        ))
    return out


# Keys that belong to the reconstruction rather than to the benchmark loop;
# used to forward the right subset when `benchmark_vs_noise` characterises.
_RECON_KEYS = frozenset({
    "p_flip_init", "n_rate_iterations", "min_confidence", "segment_blocks",
    "min_detectability", "min_contrast", "model", "emission",
    "segment_charge_jumps", "model_ramp_resets", "ramp_period",
    "n_profile_windows", "n_segment_iterations", "decoder",
})


def implied_rate_hz(reports, reported_rate_hz: float | None = None) -> float:
    """True tunnelling rate implied by what the pipeline reported.

    The pipeline does not report the truth: it loses real flips (efficiency)
    and admits spurious ones (purity), so a device at true rate ``r`` reports
    ``r * efficiency(r) / purity(r)``. A :func:`sweep_rate` measures both of
    those as functions of ``r``, so the observed count can be inverted against
    the sweep -- and the answer is then consistent with that very curve.

    This is what used to be drawn on the figure as a "measured rate" marker.
    It was removed because one trace admits three different rates and a single
    labelled line could not say which was meant:

    * ``fidelity.rate_hz`` -- posterior-threshold crossings, counted before any
      ``min_confidence`` filter, and from a different rule than the Viterbi
      flip list. Measured on a 1 s trace at contrast 2.4: 22 Hz.
    * the reported count, ``n_flips / duration`` -- same decoder and same cut
      as the sweep, but depressed by every flip the pipeline missed: 16 Hz.
    * what those imply about the device, which is this function and the only
      one on the same footing as the x-axis: 17.1 Hz.

    Parameters
    ----------
    reports : list of BenchmarkReport
        Output of :func:`sweep_rate`.
    reported_rate_hz : float, optional
        The rate to invert. Defaults to the measured trace's own reported
        count, ``n_flips / duration``.

    Returns
    -------
    float, or nan when the reported rate falls outside the swept range, where
    there is nothing to invert against and extrapolating would be a guess.
    """
    reports = list(reports)
    if not reports:
        raise ValueError("no reports to invert")
    if reported_rate_hz is None:
        f = reports[0].fidelity
        reported_rate_hz = (f.n_flips / f.duration if f.duration > 0
                            else float("nan"))
    rate = np.array([np.mean([t.rate_injected for t in r.trials])
                     for r in reports], dtype=float)
    eff = np.array([r.efficiency_clustered[0] for r in reports], dtype=float)
    pur = np.array([r.purity_clustered[0] for r in reports], dtype=float)
    ok = np.isfinite(rate) & np.isfinite(eff) & np.isfinite(pur) & (pur > 0)
    if ok.sum() < 2 or not np.isfinite(reported_rate_hz) or reported_rate_hz <= 0:
        return float("nan")
    rate, eff, pur = rate[ok], eff[ok], pur[ok]
    order = np.argsort(rate)
    rate = rate[order]
    predicted = rate * eff[order] / pur[order]
    if not np.all(np.diff(predicted) > 0):
        keep = np.concatenate([[True], np.diff(predicted) > 0])
        predicted, rate = predicted[keep], rate[keep]
    if predicted.size < 2 or not (predicted[0] <= reported_rate_hz <= predicted[-1]):
        return float("nan")
    return float(np.exp(np.interp(np.log(reported_rate_hz), np.log(predicted),
                                  np.log(rate))))


def burst_n50(points) -> float:
    """Burst multiplicity at which detection efficiency crosses one half.

    The single summary number of a :func:`sweep_burst_size` curve. Returned
    rather than annotated on the figure. ``nan`` if the curve does not cross.
    """
    pts = sorted(points, key=lambda p: p.n_qp_expected)
    q = np.array([p.n_qp_expected for p in pts], dtype=float)
    e = np.array([p.efficiency for p in pts], dtype=float)
    cross = np.flatnonzero((e[:-1] < 0.5) & (e[1:] >= 0.5))
    if not cross.size:
        return float("nan")
    i = int(cross[0])
    span = e[i + 1] - e[i]
    return float(q[i] + (0.5 - e[i]) * (q[i + 1] - q[i]) / span if span
                 else q[i])
