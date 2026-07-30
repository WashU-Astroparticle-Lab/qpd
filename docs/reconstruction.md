# Reconstructing charge-parity dynamics from a readout trace

**Methodology notes for `qpd.reconstruction`.**

These notes are written to be readable by someone who has taken an
undergraduate course in probability and linear algebra and has never seen a
hidden Markov model before. Every symbol is defined where it first appears, and
the statistical machinery is built up from scratch in [§5](#5-the-hidden-markov-model-hmm).

> **A note on the acronym.** The core of this method is a **HMM** — a *Hidden
> Markov Model*. That is a deterministic, exact inference algorithm for
> sequences. It is **not** HMC (*Hamiltonian Monte Carlo*), which is a random
> sampling method used in Bayesian fitting. They are unrelated; nothing here
> samples anything.

---

## Contents

1. [The physical problem](#1-the-physical-problem)
2. [Notation](#2-notation)
3. [Why the parity is visible at all](#3-why-the-parity-is-visible-at-all)
4. [Reducing the plane to a line](#4-reducing-the-plane-to-a-line)
5. [The hidden Markov model (HMM)](#5-the-hidden-markov-model-hmm)
6. [Case A — constant offset charge](#6-case-a--constant-offset-charge)
7. [Case B — swept offset charge](#7-case-b--swept-offset-charge)
8. [Nuisance 1 — sawtooth ramp resets](#8-nuisance-1--sawtooth-ramp-resets)
9. [Nuisance 2 — offset-charge jumps](#9-nuisance-2--offset-charge-jumps)
10. [Estimating the tunnelling rate](#10-estimating-the-tunnelling-rate)
11. [From a decoded sequence to flip times](#11-from-a-decoded-sequence-to-flip-times)
12. [Scoring against truth](#12-scoring-against-truth)
13. [What sets the limits](#13-what-sets-the-limits)
14. [Reproducing the numbers](#14-reproducing-the-numbers)
15. [Glossary](#15-glossary)

---

## 1. The physical problem

A QPD transmon is coupled to a readout resonator. The transmon's **charge
parity** — whether the number of electrons on the island is even or odd —
shifts the resonator frequency by a small amount. Because we drive the
resonator at a fixed frequency and record the transmitted microwave signal, that
frequency shift shows up as a shift in the complex transmission $S_{21}$.

The parity is not static. Every time a **quasiparticle** (a broken Cooper pair)
tunnels across the junction, the parity flips: even → odd → even → … The parity
is therefore a *random telegraph signal*: it holds a value for a random dwell
time, then toggles.

**The measurement.** We record the in-phase and quadrature voltages $I(t)$,
$Q(t)$ at a uniform sampling rate (100 kHz throughout, so one sample every
$\Delta t = 10\ \mu\text{s}$). We combine them into one complex number per
sample, $z_k = I_k + \mathrm{i}Q_k$.

**The goal.** Recover the **time of every tunnelling event** — a list of flip
times $\{\hat{t}_1, \hat{t}_2, \dots\}$ — from the trace alone.

**"Blind."** Throughout, the reconstruction is given *only* the trace and the
sampling rate. It is not told $E_J$, $E_C$, the coupling $g$, the resonator
parameters, the noise level, or how the offset charge was driven. Everything it
needs is estimated from the data. This is deliberate: it is what a reconstruction
applied to real measured data must do, and it prevents us from accidentally
scoring a method that only works because we handed it the answer.

Why this is hard, in one number: at the reference operating point a *single
sample* separates the two parity states by only about **2.7 standard
deviations of noise** — and for part of the time by much less. Classifying each
sample independently is hopeless. The whole method is about combining
information across many samples in the statistically correct way.

---

## 2. Notation

| Symbol | Meaning | Units |
|---|---|---|
| $t$ | time | s |
| $\Delta t$ | sampling period ($10\ \mu$s here) | s |
| $n$ | number of samples in the trace | — |
| $k$ | sample index, $k = 1 \dots n$; $t_k = (k-1)\Delta t$ | — |
| $z_k = I_k + \mathrm{i}Q_k$ | complex readout sample | a.u. |
| $n_g$ | offset charge, in units of one Cooper pair ($2e$) | — |
| $S_k \in \{A, B\}$ | hidden parity branch at sample $k$ | — |
| $\chi_{\rm even}, \chi_{\rm odd}$ | dispersive shift for each parity | Hz |
| $\Gamma$ | quasiparticle tunnelling rate (per state) | Hz |
| $p = \Gamma\,\Delta t$ | probability of a flip between consecutive samples | — |
| $\sigma$ | noise standard deviation, per quadrature | a.u. |
| $x_k$ | the trace projected onto one real axis (§4) | a.u. |
| $c(t)$ | *common mode*: midpoint of the two branches | a.u. |
| $h(t)$ | *signed splitting*: separation of the two branches | a.u. |
| $C = \lvert h\rvert/\sigma$ | **contrast** — separation in units of noise | — |
| $P$ | *fold period* — see §7 | s |
| $T$ | sawtooth ramp period — see §8 | s |
| $\delta$ | size of an offset-charge jump, in Cooper pairs | — |

Two conventions used everywhere:

- **Frequencies and rates are in Hz**, times in seconds (the repo-wide convention).
- The Cooper-pair box is **periodic in $n_g$ with period 1**. One Cooper pair
  is $2e$, so "half a Cooper pair" means $\Delta n_g = 0.5$.

---

## 3. Why the parity is visible at all

Two facts about the physics do all the work.

**Fact 1 — odd parity is even parity, shifted by half a Cooper pair.**

$$\chi_{\rm odd}(n_g) \;=\; \chi_{\rm even}\!\left(n_g + \tfrac{1}{2}\right)$$

and $\chi_{\rm even}$ is periodic with period 1 and is an even function of $n_g$.
So if we define the **signed splitting**

$$h(n_g) \;\equiv\; \chi_{\rm odd}(n_g) - \chi_{\rm even}(n_g),$$

then $h$ is periodic and, crucially, it **passes through zero**. The zeros sit at

$$n_g = 0.25 \pmod{0.5}.$$

At those offset charges the two parity states produce *exactly the same*
resonator frequency. They are called **parity-blind points**: no measurement,
however long, can distinguish the parities there. This is physics, not a
shortcoming of the algorithm, and it recurs throughout these notes.

Numerically, for the reference device ($E_J/E_C = 12$, $g = 150$ MHz):

| $n_g$ | $\chi_{\rm even}$ (MHz) | $\chi_{\rm odd}$ (MHz) | $h$ (MHz) |
|---|---|---|---|
| 0.00 | 27.94 | 22.48 | −5.46 |
| 0.10 | 27.28 | 22.89 | −4.39 |
| **0.25** | **24.85** | **24.85** | **0.00** |
| 0.40 | 22.89 | 27.28 | +4.39 |
| 0.50 | 22.48 | 27.94 | +5.46 |

Note the sign flip either side of $n_g = 0.25$. Which branch is "upper" swaps
every time we cross a blind point — a bookkeeping detail that becomes important
in §7.

**Fact 2 — the splitting is small compared with the resonator linewidth.**

The two branches shift the resonator by $\lesssim 5.5$ MHz, against a linewidth
$\kappa \approx 243$ kHz but with the drive parked far out on the Lorentzian
tail. The practical consequence is that as the parity flips, the measured
$S_{21}$ moves back and forth along an essentially **straight line segment** in
the complex plane. It does not swing around an arc. (Measured on the reference
device, the direction of the difference vector is constant to about 3 mrad
across the whole sweep.)

Fact 2 is what makes §4 possible.

---

## 4. Reducing the plane to a line

Each sample $z_k$ is a point in the complex plane. Under branch $S_k$ its
expected position is $\mu_{S_k}(t_k)$, and the measurement adds isotropic
Gaussian noise of width $\sigma$ **independently to $I$ and to $Q$**:

$$z_k \;=\; \mu_{S_k}(t_k) \;+\; \varepsilon_k, \qquad
\varepsilon_k \sim \mathcal{N}(0,\sigma^2) + \mathrm{i}\,\mathcal{N}(0,\sigma^2).$$

By Fact 2 the two means $\mu_A$ and $\mu_B$ differ only along one fixed
direction. Call the unit vector along that direction $\hat{d}$ (a complex number
of modulus 1). Then all the parity information lies in the component of $z_k$
along $\hat{d}$; the perpendicular component is pure noise and can be discarded.

Define the **projection**

$$x_k \;\equiv\; \operatorname{Re}\!\left[\, \overline{\hat{d}}\,(z_k - o) \,\right],$$

where $\overline{\hat{d}}$ is the complex conjugate and $o$ is any convenient
origin (we use the mean of the trace). Because a projection of an isotropic
Gaussian is a 1-D Gaussian of the *same* width, $x_k$ satisfies

$$x_k \sim \mathcal{N}\!\big(m_{S_k}(t_k),\ \sigma^2\big),$$

with $m_A, m_B$ the projected branch means. **The two-dimensional problem is now
one-dimensional, with no loss of information.**

**Finding $\hat{d}$ blind.** The cloud of all samples is an isotropic blob
smeared out along one line. Its covariance matrix therefore has a large
eigenvalue along $\hat{d}$ and a small one perpendicular to it. So:

- $\hat{d}$ = principal (major) eigenvector of the sample covariance;
- $\sigma$ = square root of the **minor** eigenvalue.

That second point is a small gift: the minor axis carries noise *only*, so it
gives an estimate of $\sigma$ that needs no assumption about the flip rate or the
ramp. (`estimate_direction` in [`emission.py`](../src/qpd/reconstruction/emission.py).)

It is convenient to write the two projected means as a midpoint plus a splitting:

$$m_A = c + \tfrac{h}{2}, \qquad m_B = c - \tfrac{h}{2},$$

so $c$ is the **common mode** and $h$ the **signed splitting**. The **contrast**

$$\boxed{\;C \;=\; \frac{|h|}{\sigma}\;}$$

is the single most important number in these notes: it is the branch separation
measured in units of the noise. At the reference operating point the median
contrast is $C \approx 2.7$.

---

## 5. The hidden Markov model (HMM)

This section assumes no prior exposure. We build the model, then derive the two
algorithms we need.

### 5.1 What we know and what we want

We know the sequence of measurements $x_1, \dots, x_n$. We want the sequence of
parities $S_1, \dots, S_n$, which we cannot see — hence *hidden*.

A single sample tells us very little: with $C \approx 2.7$ the two candidate
Gaussians overlap substantially, so a per-sample guess is wrong maybe 10% of the
time — which over $10^6$ samples would mean $10^5$ spurious "flips." What saves
us is that **parity flips are rare**. At $\Gamma = 10$ Hz the parity holds for
about 100 ms, i.e. $10^4$ consecutive samples. Combining $10^4$ samples that each
carry contrast 2.7 gives an overwhelming determination of which branch we are on.

The HMM is exactly the machinery for combining evidence across a sequence in a
way that respects "flips are rare."

### 5.2 The two ingredients

**(a) The transition model** (the prior on how the hidden state evolves).

The parity is a *Markov chain*: the probability of the next state depends only on
the current state, not on the whole history. For a symmetric telegraph with rate
$\Gamma$ per state, the probability of flipping between two consecutive samples is

$$p \;=\; \Gamma\,\Delta t .$$

(At $\Gamma = 10$ Hz and $\Delta t = 10\ \mu$s, $p = 10^{-4}$ — flips are indeed
rare per sample.) The transition matrix is

$$\mathbf{T} \;=\; \begin{pmatrix} 1-p & p \\ p & 1-p \end{pmatrix},
\qquad T_{ss'} = P(S_{k+1}=s' \mid S_k = s).$$

This is the term that stops one noisy sample from being read as a flip: flipping
"costs" a factor $p/(1-p)$ in probability, so the evidence has to be worth more
than that.

**(b) The emission model** (how a hidden state produces an observation).

Given the branch, the sample is Gaussian about that branch's mean:

$$P(x_k \mid S_k = s) \;=\; \frac{1}{\sqrt{2\pi\sigma^2}}
\exp\!\left[-\frac{\big(x_k - m_s(t_k)\big)^2}{2\sigma^2}\right].$$

We work with the logarithm, $\ell_k(s) \equiv \log P(x_k \mid S_k=s)$, because
products of $10^6$ small numbers underflow any floating-point format.

Note the means carry a time argument: in the ramped case (§7) they *move* from
sample to sample. The HMM does not care — it only ever needs $\ell_k(s)$.

### 5.3 The forward–backward algorithm (posterior probabilities)

We want, for each $k$, the probability that the branch was $B$ **given the whole
trace**:

$$\gamma_k \;\equiv\; P\!\left(S_k = B \,\middle|\, x_1,\dots,x_n\right).$$

Note "given the whole trace" — including samples *after* $k$. Data from the
future genuinely helps: if the next 5000 samples all look like branch $B$, that
is strong evidence sample $k$ was already $B$.

Define two quantities:

$$\alpha_k(s) \;\equiv\; P\!\left(x_1,\dots,x_k,\ S_k=s\right)
\qquad\text{(forward)}$$
$$\beta_k(s) \;\equiv\; P\!\left(x_{k+1},\dots,x_n \,\middle|\, S_k=s\right)
\qquad\text{(backward)}$$

**Forward recursion.** To be in state $s'$ at step $k+1$ having seen the data so
far, we must have been in *some* state $s$ at step $k$, transitioned, and then
emitted $x_{k+1}$:

$$\boxed{\;\alpha_{k+1}(s') \;=\; \Big[\textstyle\sum_{s} \alpha_k(s)\,T_{ss'}\Big]\; P(x_{k+1}\mid s')\;}$$

started from $\alpha_1(s) = P(S_1=s)P(x_1\mid s)$, with $P(S_1=s) = 1/2$ (we have
no prior preference for even or odd).

**Backward recursion.** Symmetrically, running from the end of the trace inwards:

$$\boxed{\;\beta_k(s) \;=\; \textstyle\sum_{s'} T_{ss'}\;P(x_{k+1}\mid s')\;\beta_{k+1}(s')\;}$$

started from $\beta_n(s) = 1$.

**Combining.** Because the chain is Markov, past and future are conditionally
independent given the present state, so

$$P(S_k = s,\ \text{all data}) \;=\; \alpha_k(s)\,\beta_k(s),$$

and normalising over the two states gives what we wanted:

$$\gamma_k \;=\; \frac{\alpha_k(B)\beta_k(B)}{\alpha_k(A)\beta_k(A) + \alpha_k(B)\beta_k(B)}.$$

Both recursions visit each sample once, so the cost is $O(n)$ — linear in the
length of the trace, with only two states to track.

**Numerical care.** $\alpha_k$ shrinks geometrically and would underflow within a
few hundred samples. The implementation rescales $\alpha$ to sum to 1 at every
step and accumulates the discarded factors as a running log. This is exact (the
scale factors cancel in $\gamma_k$) and the accumulated log is precisely the
total log-likelihood $\log P(x_1,\dots,x_n)$, which we reuse in §8 for model
selection. See `forward_backward` in [`hmm.py`](../src/qpd/reconstruction/hmm.py).

**How to read $\gamma_k$.** It is a soft curve between 0 and 1. At a real flip it
swings sharply from ~0 to ~1 over a few samples. Inside a parity-blind window,
where the data says nothing, it drifts back toward $1/2$ — the model reports
honest ignorance instead of guessing.

### 5.4 The Viterbi algorithm (the single best sequence)

The posterior gives a per-sample probability. To *count events* we want one
definite sequence: the single most likely assignment $S_1,\dots,S_n$ as a whole.
That is the **Viterbi** algorithm — dynamic programming, structurally the same as
the forward recursion but with $\sum$ replaced by $\max$:

$$\boxed{\;\delta_{k+1}(s') \;=\; \max_{s}\Big[\delta_k(s) + \log T_{ss'}\Big] \;+\; \ell_{k+1}(s')\;}$$

working in logs throughout ($\delta_k(s)$ is the log-probability of the best path
that ends in state $s$ at step $k$). At each step we record *which* $s$ achieved
the maximum; at the end we take the better final state and walk the recorded
choices backwards to read off the whole path. Also $O(n)$.

The result is a clean square wave. Flip events are exactly its transitions.

### 5.5 Why this beats a threshold

A naive detector — smooth the trace, threshold it, call each crossing a flip —
has a fixed averaging window, and you must choose it. Too short and noise
manufactures flips; too long and real flips are smeared away.

The HMM has no such knob. Weighting the evidence by the dwell-time prior makes it
behave like a **matched filter whose averaging length adapts automatically**: it
integrates for a long time where the branches are close (low contrast) and
responds quickly where they are well separated. That adaptivity is what makes a
per-sample contrast of order unity workable at all.

---

## 6. Case A — constant offset charge

Start with the easy regime, because it isolates the HMM from everything else.

If $n_g$ is **held fixed**, the two branch means do not move. In the I/Q plane
you see two stationary blobs, and the trace hops between them. There is no fold
period, no blind point sweeping past, no ramp reset. Entry point:
`reconstruct_parity_flips_static` ([`static_bias.py`](../src/qpd/reconstruction/static_bias.py)).

The projected data is a **two-component Gaussian mixture**:

$$P(x) \;=\; w\,\mathcal{N}(x; m_A, \sigma^2) \;+\; (1-w)\,\mathcal{N}(x; m_B, \sigma^2),$$

with $w$ the fraction of time spent in branch $A$. We fit $m_A, m_B, \sigma, w$
by **expectation–maximisation (EM)**, which alternates two steps until converged:

- **E-step** — given current parameters, compute each sample's *responsibility*,
  the probability it came from branch $A$:
  $$r_k \;=\; \frac{w\,\mathcal{N}(x_k; m_A,\sigma^2)}{w\,\mathcal{N}(x_k; m_A,\sigma^2) + (1-w)\,\mathcal{N}(x_k; m_B,\sigma^2)}.$$
- **M-step** — re-estimate the parameters as responsibility-weighted averages:
  $$m_A = \frac{\sum_k r_k x_k}{\sum_k r_k},\qquad
    m_B = \frac{\sum_k (1-r_k) x_k}{\sum_k (1-r_k)},\qquad
    w = \frac{1}{n}\sum_k r_k,$$
  and $\sigma^2$ as the pooled weighted variance.

EM is initialised by splitting at the **median** of $x$. That is a deliberate
choice: splitting at the extremes (k-means style) biases the two centres outward
when the blobs overlap, whereas for the balanced occupancy of a symmetric
telegraph the median split is unbiased.

Then feed $m_A$, $m_B$, $\sigma$ into the HMM of §5 and read off the transitions.
That is the entire static pipeline.

### 6.1 A trap: EM always "succeeds"

Fit two Gaussians to data that is really *one* Gaussian and EM will still return
two components — it splits the single blob into a spurious pair, typically
reporting a contrast near 0.9. **So the fitted contrast alone cannot tell you
whether the bias point carries any parity information.** At the parity-blind
charge $n_g = 0.25$, where the truth is that no information exists, the fit looks
superficially the same as at a good bias point.

The honest tell is the *combination* of contrast and decoded dwell time. Define

$$\text{detectability} \;=\; C\,\sqrt{N_{\rm dwell}},
\qquad N_{\rm dwell} = 1/p,$$

the contrast integrated over one dwell (the $\sqrt{N}$ is the usual averaging
gain for $N$ independent samples). If the decoder is really segmenting noise, the
fitted rate is inflated, the dwell collapses, and the product falls. On the
reference device every working bias point scores $\geq 100$ and every failing one
$\leq 47$, with nothing in between — so the threshold at 70 sits in an empty gap.
Traces below it are flagged `degenerate`, and the correct response is to **move
the bias point**, not to salvage the output.

---

## 7. Case B — swept offset charge

Now the case the notebook actually simulates: $n_g$ is **ramped** with a sawtooth
at 500 Hz, spanning 5.5 Cooper pairs. Entry point
`reconstruct_parity_flips_ramped` ([`reconstruct.py`](../src/qpd/reconstruction/reconstruct.py)).

Sweeping helps and hurts. It helps because no flip is stuck at a bad bias point
forever — the contrast sweeps, so every flip eventually sees good separation.
It hurts because the branch means now *move*, and we must learn how.

### 7.1 The fold period

Since $|h|$ is periodic in $n_g$ with period 0.5, and a linear ramp makes $n_g$
proportional to $t$, the splitting is periodic **in time**:

$$P \;=\; \frac{0.5}{\text{ramp slope}} \qquad\text{(the \emph{fold period})}.$$

For the reference scenario the slope is $5.5 \times 500 = 2750$ per second, so
$P = 181.8\ \mu$s — about 18 samples.

### 7.2 Finding $P$ blind: look at the variance, not the mean

Here is the key trick. The **common mode barely moves** — the midpoint $c$
between the branches is nearly constant. What the ramp modulates is the
*separation*. So the parity signal lives in the trace's **spread**, not its mean,
and averaging the trace would destroy it.

Consider the square of the projection. Over a phase bin where the parity is
equally likely to be $A$ or $B$, the mean of $x$ is $c$ and

$$\mathbb{E}\!\left[(x-c)^2\right]
\;=\; \underbrace{\sigma^2}_{\text{noise}} \;+\; \underbrace{\left(\tfrac{h}{2}\right)^2}_{\text{branch offset}}
\;=\; \sigma^2 + \frac{h^2}{4}.$$

*Derivation:* with probability $1/2$ the sample sits at $c + h/2$ plus noise, and
with probability $1/2$ at $c - h/2$ plus noise. In either case the offset from
$c$ is $\pm h/2$, whose square is $h^2/4$; the noise contributes $\sigma^2$ and is
independent, so the variances add.

So $x^2$ contains a strong periodic tone at the fold frequency $1/P$. We locate it
with a periodogram (an FFT of $x^2$), then **refine** it, because raw FFT
resolution is nowhere near good enough: a 5 s trace at $P = 182\ \mu$s contains
27 500 cycles, so a relative period error of only $3\times10^{-5}$ slides the
model three quarters of a cycle out of phase by the end of the trace. Refinement
maximises the *folded contrast* — fold the trace at a trial period, bin by phase,
and take the variance of the resulting profile; this is sharpest at the true
period because a wrong period smears the profile flat.

### 7.3 Folding, and the profile

With $P$ in hand, assign each sample a **phase**
$\varphi_k = \operatorname{frac}(t_k/P) \in [0,1)$, bin the samples by $\varphi$,
and in each bin compute the mean and variance of $x$. Then

$$c(\varphi) = \text{bin mean},\qquad
|h(\varphi)| = 2\sqrt{\max\big(\text{bin variance} - \sigma^2,\ 0\big)},$$

inverting the relation above. Because every bin pools thousands of samples from
across the whole trace, the profile is measured very precisely even though any
individual sample is noisy.

### 7.4 Recovering the sign

The fold gives $|h|$, but the HMM needs the *signed* $h$ — it must know which
branch is currently upper. From §3, $h$ changes sign at every parity-blind point,
and there is exactly one blind point per fold period (the minimum of $|h|$). So
the sign simply **alternates from one half-cycle to the next**:

$$h(t) \;=\; (-1)^{\lfloor t/P - \varphi_{\rm blind}\rfloor}\,\big|h(\varphi)\big|.$$

The overall sign is not identifiable — flipping it globally just swaps the labels
"even" and "odd." That is harmless: it does not move any flip *time*, which is
all we are reporting.

### 7.5 Robustness of the period estimate

One refinement worth knowing, because it caused a real failure during
development. If a charge jump (§9) sits in the middle of the trace, the
whole-trace folding objective can trade the phase step against a small period
error, biasing the period by $3\times10^{-5}$ — enough to ruin everything
downstream. The fix: estimate the period **independently in several windows and
take the median**. A jump falls in only one window, so the median is untouched by
it. (`estimate_fold_period`, `min_cycles_per_window`.)

---

## 8. Nuisance 1 — sawtooth ramp resets

This is the dominant systematic in the swept case, and the most instructive.

The sawtooth spans 5.5 Cooper pairs. In units of half-pairs that is
$5.5/0.5 = \mathbf{11}$ — an **odd** number. So when the ramp resets, $n_g$ jumps
by 5.5, which modulo 1 is exactly **0.5**: half a Cooper pair. And by Fact 1
(§3), shifting $n_g$ by half a pair maps each parity branch *exactly onto the
other*.

**A ramp reset is therefore observationally identical to a parity flip.** The
signal does precisely what a tunnelling event would do. No per-event test can
distinguish them, because there is no difference to detect.

The numbers make this severe. At a 500 Hz ramp there are 500 resets per second;
at $\Gamma = 10$ Hz there are ~10 real flips per second. That is **~50 artefacts
for every real event**.

### 8.1 Why you cannot fix this afterwards

The tempting approach is to detect all events, notice that some are strictly
periodic, and delete those. This fails badly here. To delete one event per 2 ms
cycle you need a matching window; with a window comparable to the timing
tolerance you would also delete a substantial fraction of the *real* flips that
happen to fall near a reset phase. You cannot win: the artefacts are 50× more
numerous, so even a small collateral rate destroys the result.

### 8.2 The fix: put the resets in the model

Resets are strictly **periodic**, while tunnelling is **Poisson** (random). That
is the one statistical difference available, and it is a property of the
*ensemble*, not of any single event.

So: run a short first-pass decode. Its event list is dominated by reset
artefacts. Fold those times at candidate periods and look for a sharp pile-up at
one phase — random events spread out, a comb piles up. The search is restricted
to integer multiples of the already-known fold period $P$ (the ramp spans a whole
number of half-pairs, so the reset interval must be $T = mP$ for integer $m$),
which makes it both cheap and sharp. On the reference scenario this returns
$m = 11$, $T = 2.0000$ ms with a pile-up significance of several hundred sigma.

Then, instead of deleting anything, **add an extra sign flip at each reset time**
to the sign schedule of §7.4. The emission model now knows that the branches swap
at those instants, so the decoder simply follows them through — and reports **no
event at all**. The artefacts never enter the output.

The comb is accepted only if it **improves the total log-likelihood** (available
free from the forward pass, §5.3). This gate means a trace with no such ramp —
a static bias, an even-span sawtooth, a triangle wave with no reset — is left
untouched.

Measured effect: hard $F_1$ goes from **0.20 → 0.94**, and **0 of 2500** resets
in a 5 s trace leak into the output.

---

## 9. Nuisance 2 — offset-charge jumps

Occasionally (~0.1 Hz) a stray charge shifts $n_g$ abruptly by some $\delta$.

A jump and a flip look different, and the difference is what lets us handle them:

- a **parity flip** moves the signal *between* two branch curves that themselves
  stay put — a transient change;
- a **charge jump** slides *both* curves along the ramp phase, by $2\delta$ in
  fold-period units (one fold period is half a Cooper pair, hence the factor 2).
  The learned model is then wrong from the jump onward — a persistent change.

So we detect the jump as a **change point** in the fitted ramp phase, re-fit one
phase offset per segment, and emit no flip at the boundary. Jump *amplitudes* are
not reported; they are not a target here.

**This cannot be skipped.** A jump near $\delta = 0.25$ shifts the phase by
exactly half a cycle, which inverts which branch is which; the decoder then
flip-flops constantly and $F_1$ collapses from 0.98 to **0.21**.

Two subtleties, both found the hard way:

1. **Chicken-and-egg.** You need a splitting profile to find the jump, but a
   whole-trace profile is already corrupted *by* the jump. Worst at
   $\delta = 0.25$, which blends the profile with its own half-period shift,
   making it nearly period-$\tfrac12$ — at which point an offset of 0 is
   indistinguishable from 0.5. Fix: build the profile from a **single window**
   (jump-free for every window the jump misses) and pick the window with the
   highest contrast, since blending always costs contrast.
2. **Coupling with the period.** A phase step biases the period (§7.5) and a
   biased period drifts the phase in a way that mimics a step. The two are
   estimated by **alternating** between them until stable.

**Irreducible case.** A jump of exactly $\delta = \pm 0.5$ (mod 1) maps each
branch onto the other — the same degeneracy as a ramp reset, but occurring at a
random time, so the periodicity argument of §8 does not apply. It is reported as
a flip. This is a genuine limit of the measurement, not of the algorithm.

---

## 10. Estimating the tunnelling rate

The transition probability $p$ is itself unknown. It is estimated by **hard EM**:

1. decode with the current $p$;
2. count transitions in the decoded path, set $p \leftarrow (\text{count})/(n-1)$;
3. repeat until stable (a handful of iterations).

The rate only sets the decoder's *prior*, so this crude scheme is adequate — a
full Baum–Welch update would cost an extra sweep per iteration for no measurable
gain. Measured recovery: 1 Hz → 1.0, 10 Hz → 9.2, 100 Hz → 98.8, 300 Hz → 292.7.

---

## 11. From a decoded sequence to flip times

The Viterbi path changes state between two samples, so the flip happened
somewhere in that $10\ \mu$s interval. Rather than taking the midpoint, we read
off where the **posterior odds cross even** — linearly interpolating $\gamma_k$
between the bracketing samples to find where it passes $1/2$. With a $10\ \mu$s
sample period against a $0.5$ ms matching tolerance this is a small correction,
but it is free and unbiased.

Each flip also carries a **confidence**: the size of the posterior swing across
the transition, in $[0,1]$. This collapses toward zero for flips inside a
parity-blind window, giving a principled handle for trading recall against
precision. (Note: in the measured failure modes the false positives turned out to
be *high*-confidence, so this is a diagnostic, not a cure-all.)

---

## 12. Scoring against truth

Reusing `qpd.mlebench.grade`, so reconstruction is scored the same way the
benchmark datasets are.

**Matching.** Predicted times are matched one-to-one to true times, minimising
total timing error subject to a hard tolerance of 0.5 ms. This is a bipartite
assignment problem, solved exactly with the Hungarian algorithm
(`scipy.optimize.linear_sum_assignment`). To stay tractable on dense trains the
tolerance-gated graph is split into connected components and solved per
component — exact, since no edge crosses a component boundary.

**Metrics.** With $N_{\rm true}$ true events, $N_{\rm pred}$ predicted, and
$N_{\rm match}$ matched:

$$\text{efficiency (recall)} = \frac{N_{\rm match}}{N_{\rm true}},\qquad
\text{purity (precision)} = \frac{N_{\rm match}}{N_{\rm pred}},$$

$$F_1 \;=\; \frac{2\,N_{\rm match}}{N_{\rm true} + N_{\rm pred}}
\;=\; \text{harmonic mean of efficiency and purity}.$$

$F_1$ is used because it punishes both failure modes: missing real events *and*
inventing fake ones. A detector that reports nothing has perfect purity; one that
reports an event every sample has perfect efficiency. Only $F_1$ rejects both.

A **soft** $F_1$ additionally weights each match by timing accuracy,
$q = \exp[-(\Delta t/\text{scale})^2]$ with scale = 0.25 ms, so a sloppily-timed
match earns partial credit.

Separately we report the **timing bias** (mean signed residual) and **rms**.
Keeping these apart from $F_1$ matters: they answer different questions —
"did we find the events?" versus "did we time them correctly?" — and in the burst
regime they behave very differently (§13).

---

## 13. What sets the limits

**(a) Contrast.** Everything is governed by $C = |h|/\sigma$. Measured on the
reference device, sweeping the noise alone:

| $\sigma$ | contrast | efficiency | purity | $F_1$ | timing rms |
|---|---|---|---|---|---|
| 0.5e−4 | 10.7 | 1.000 | 1.000 | 1.000 | 4 µs |
| 2e−4 (ref) | 2.7 | 1.000 | 1.000 | 1.000 | 16 µs |
| 3e−4 | 1.8 | 0.993 | 0.979 | 0.986 | 27 µs |
| 5e−4 | 1.1 | 0.985 | 0.898 | 0.939 | 51 µs |
| 8e−4 | 0.7 | 0.601 | 0.011 | 0.021 | **calibration fails** |

Degradation is graceful down to $C \approx 1$; timing rms scales roughly as
$1/C$, flooring near the sample period. The collapse at $C \approx 0.7$ is a
**calibration** failure, not a decoding one: the reset comb of §8 no longer
stands out of the first-pass event list, and the 500 Hz artefacts flood the
output. Useful design rule: rms $\approx (\Delta t/2)\times(5/C)$, so 25 µs timing
needs $C \gtrsim 2$.

**(b) Event crowding.** Inside a quasiparticle burst, tens of tunnels arrive
within milliseconds. Two flips separated by less than about two samples (20 µs)
cannot be resolved at 100 kHz *by anything* — the parity returns to its original
value before the next sample. This sets a ceiling:

| flips in window | efficiency | timing bias | unresolvable pairs |
|---|---|---|---|
| 1 | 1.000 | −0.6 µs | 0.000 |
| 7–10 | 0.889 | −1.5 µs | 0.028 |
| 17–25 | 0.794 | −1.3 µs | 0.070 |
| 41+ | 0.561 | −1.6 µs | 0.141 |

Two things to read here. First, **timing bias stays at −1 to −2 µs throughout**:
what crowding costs you is *missed events*, not mistimed ones. Second, at the
highest crowding the achieved efficiency (0.56) is well below the intrinsic
ceiling (~0.86), so part of that loss *is* the algorithm: the globally-fitted
dwell prior (10–40 Hz) is far too slow for in-burst rates near 3 kHz, so the HMM
under-segments. A burst-aware HMM — hidden state = (branch) × (quiet/burst
regime), with a much larger $p$ in the burst regime — would close the gap.

**(c) Parity-blind points.** Flips landing where $h \approx 0$ are unrecoverable
in principle. In the swept case these windows are brief and frequent, so little
is lost; in the static case a bad bias point loses everything (§6.1).

**(d) The $\delta = \pm 0.5$ degeneracy.** §9. Irreducible in a single trace.

---

## 14. Reproducing the numbers

```bash
# regression gate (~1 min)
python checks/check_parity_reconstruction.py

# performance studies: rate / noise / burst-crowding / device sweeps
python checks/study_parity_reconstruction.py            # all
python checks/study_parity_reconstruction.py noise      # one
```

A worked end-to-end demonstration, with the learned branch means overlaid on the
trace, the posterior, and the reconstruction-vs-truth comparison, is §7 of
[`notebooks/simulation.ipynb`](../notebooks/simulation.ipynb).

Minimal use:

```python
from qpd.reconstruction import reconstruct_parity_flips_ramped, score_flips

rec = reconstruct_parity_flips_ramped(result.iq, sample_rate=1e5)
rec.flip_times          # reconstructed tunnelling times [s]
rec.reset_comb          # the ramp-reset schedule it found
rec.charge_jump_times   # jumps found and corrected as nuisances

score_flips(result.flip_times, rec.flip_times)   # efficiency, purity, F1, rms
```

---

## 15. Glossary

**Blind (reconstruction).** Using only the recorded trace and its sample rate;
no device, resonator, or noise parameters supplied.

**Contrast $C$.** Branch separation in units of the noise standard deviation,
$|h|/\sigma$. The controlling figure of merit.

**Common mode $c$.** The midpoint between the two branch means. Nearly constant,
which is why the parity signal lives in the variance, not the mean.

**Dwell time.** How long the parity holds before the next tunnelling event;
mean $1/\Gamma$.

**Efficiency / purity.** Recall and precision for event detection: fraction of
true events found, and fraction of reported events that are real.

**EM (expectation–maximisation).** Iterative fitting for models with hidden
labels: alternate between assigning soft responsibilities and re-estimating
parameters.

**Emission model.** In an HMM, the probability of an observation given the hidden
state. Here: a Gaussian about the current branch mean.

**Fold period $P$.** The period with which the branch splitting repeats in time
under a linear $n_g$ ramp; $P = 0.5/\text{slope}$.

**Forward–backward.** The exact $O(n)$ algorithm giving the posterior
probability of each hidden state given the *entire* sequence.

**HMM (hidden Markov model).** A model of a hidden state sequence with Markov
transitions, observed only through noisy emissions. *Not* HMC (Hamiltonian Monte
Carlo), which is unrelated.

**Parity-blind point.** An offset charge, $n_g = 0.25 \pmod{0.5}$, where the two
parity branches produce identical signals and no measurement can separate them.

**Projection.** Collapsing the complex I/Q sample onto the single real axis along
which the two branch means differ; lossless here, since the perpendicular
direction is pure noise.

**Ramp reset.** The sawtooth's return to the start of its range. Because the span
is an odd number of half Cooper pairs, it swaps the two branches — identical in
the data to a parity flip.

**Telegraph (random telegraph signal).** A two-state process that holds each
state for a random (exponential) time before toggling.

**Viterbi.** Dynamic programming that returns the single most likely hidden
state *sequence*, as opposed to per-sample probabilities.
