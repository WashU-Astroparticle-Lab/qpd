# From Eq. (6) to Eq. (7): the EPR zero-point-fluctuation identity

**Source:** Section II.A (transmon coupled to a cavity) of [Minev et al., arXiv:2010.00620](https://arxiv.org/abs/2010.00620).

> Equations are pre-rendered to images because Warp's markdown viewer has no math engine.
> Inline symbols use Unicode ($\varphi$, $\omega$, $\hbar$, $\hat a$, $n$).

---

## What the two equations are

**Eq. (6)** — the *quantum* definition of the participation ratio: the fraction of **inductive** energy stored in the junction when only mode $m$ is excited (state $|\psi_m\rangle$ = Fock/coherent excitation of mode $m$), with **normal ordering**:

![Eq 6](img/02_eq06.svg)

**Eq. (7)** — the result we want to reach:

![Eq 7](img/02_eq07.svg)

The goal is to show four facts collapse Eq. (6) into Eq. (7).

---

## Step 1 — Junction flux in the mode basis

The junction reduced flux is linear in the mode operators; isolating the single excited mode $m$, and writing its linear inductive energy ($L_J = \varphi_0^2/E_J \Rightarrow \Phi_J^2/2L_J = \tfrac12 E_J \hat\varphi_J^2$):

![flux and inductive energy](img/02_flux.svg)

Here **$\varphi_m \equiv \varphi_{mJ}$ is exactly the zero-point fluctuation amplitude we are solving for.**

---

## Step 2 — Numerator: junction energy, normal-ordered

Normal ordering gives $\langle n|{:}(\hat a+\hat a^\dagger)^2{:}|n\rangle = 2n$ (no vacuum +1):

![normal ordering identity](img/02_normord.svg)

so the junction's inductive-energy expectation is

![numerator](img/02_num.svg)

---

## Step 3 — Denominator: equipartition

A harmonic mode splits its energy equally between inductive (magnetic) and capacitive (electric) parts, $\langle\hat\Phi^2/2L\rangle = \langle\hat Q^2/2C\rangle$. The normal-ordered total energy of mode $m$ is $\hbar\omega_m\cdot n_m$, so the **total inductive** energy is half of that:

![denominator](img/02_den.svg)

---

## Step 4 — Take the ratio (the n_m cancels)

![ratio and result](img/02_ratio.svg)

The excitation number $n_m$ cancels — the participation is **intensive**, as a "fraction of energy" must be. Inverting gives Eq. (7).

---

## Reminder: what is normal ordering?

**Normal ordering**, written `:Ô:`, reorders a product of ladder operators so that **all creation operators $\hat a^\dagger$ sit to the left of all annihilation operators $\hat a$** — reshuffling them *as if they commuted* (i.e. ignoring $[\hat a,\hat a^\dagger]=1$ during the move). Concretely:

![normal-ordering definition](img/02_no_def.svg)

The point is the dropped commutator. For the square that appears in our derivation:

![normal-ordered square vs ordinary](img/02_no_sq.svg)

So $:(\hat a+\hat a^\dagger)^2:$ is exactly the ordinary $(\hat a+\hat a^\dagger)^2$ **minus the vacuum term** that the commutator $[\hat a,\hat a^\dagger]=1$ produces. That "+1" is the zero-point fluctuation; normal ordering throws it away.

Two ways to see what it does:

- **Fock state:** $\langle n|{:}(\hat a+\hat a^\dagger)^2{:}|n\rangle = 2n$, versus $\langle n|(\hat a+\hat a^\dagger)^2|n\rangle = 2n+1$.
- **Coherent state:** taking the expectation of a normal-ordered operator is the same as the **classical replacement** $\hat a \to \alpha$, $\hat a^\dagger \to \alpha^*$:

![normal ordering = classical replacement](img/02_no_coh.svg)

This last property is the deep reason the EPR method uses it: normal ordering is the bridge that makes the *quantum* junction-energy expectation reduce to the *classical* field energy that the finite-element solver computes — with the divergent vacuum piece removed.

---

## Is normal ordering an approximation? (No — it's exact)

A natural worry: by reshuffling operators "as if they commuted," are we ignoring $[\hat a,\hat a^\dagger]=1$ and making an approximation? **No.** Three reasons.

**(1) It's a definition, not an approximation.** `:Ô:` is a *new operator*, equal to $\hat O$ plus a known c-number — the commutator term is **separated out, not deleted**:

![operator = normal-ordered + c-number](img/02_decomp.svg)

We are only choosing *which* operator to call "the junction energy." Nothing is approximated.

**(2) For the real Josephson cosine, normal ordering is an exact rewriting.** With $\hat\varphi = \varphi_{\rm zpf}(\hat a+\hat a^\dagger)$, the exact BCH identity (valid because $[\hat a,\hat a^\dagger]$ is a c-number) gives

![BCH / displacement identity](img/02_bch.svg)

and therefore, with **no truncation**:

![exact cosine identity](img/02_exact_cos.svg)

The commutator is not ignored — it *generates* the prefactor $e^{-\varphi_{\rm zpf}^2/2} < 1$, which is the real, vacuum-fluctuation reduction of the effective Josephson energy (a Debye–Waller-like factor). It's a re-summation, not a dropped term.

**(3) What we drop is unobservable; observable vacuum physics is kept.** The set-aside piece (the +1, the $e^{-\varphi_{\rm zpf}^2/2}$ baseline) is the **zero-point** contribution. Excluding it from the *participation ratio* is legitimate because it (a) is an additive constant $\to$ no dynamics, (b) **cancels** in every measurable energy *difference*, and (c) would **diverge** if summed over all modes. Meanwhile the *nontrivial* vacuum effects are retained — the **Lamb shift** $\Delta_m$ in $\hat H_{\rm full}$ ($\omega_m' = \omega_m - \Delta_m$) is itself a zero-point effect and is kept.

> **Bottom line:** normal ordering is an exact operator identity. The commutator is accounted for —
> as the explicit c-number we subtract, or (for the cosine) as the exact $E_J$ renormalization
> $e^{-\varphi_{\rm zpf}^2/2}$. This is also why EPR predictions match experiment to a few percent: no real
> physics is thrown away.

---

## Why normal ordering is essential (the crux)

Without `: :`, one has $\langle n|(\hat a+\hat a^\dagger)^2|n\rangle = 2n+1$. The extra +1 is the vacuum contribution, and it would:

1. make $p_m$ depend on the excitation level $n_m$ (so it's no longer a clean "fraction"), and
2. sum to a **divergent** zero-point energy over infinitely many modes.

Normal ordering subtracts the vacuum, leaving the participation of the *excitation* — an intensive, $n_m$-independent number. (Paper: Supplementary Section A6.)

---

## One-line physical reading

Of the **$\tfrac12\hbar\omega_m$** of inductive energy carried by *one added photon* in mode $m$, the fraction **$p_m$** lives in the junction. That junction share equals **$E_J \varphi_m^2$**, which rearranges directly to

> **$\varphi_m^2 = p_m \cdot \hbar\omega_m / (2 E_J)$.**

This is the bridge from the *classical* eigenmode solve (which gives $p_m$) to the *quantum* zero-point fluctuations $\varphi_m$ that set every nonlinear coupling. The general multi-junction version is Eq. (21).
