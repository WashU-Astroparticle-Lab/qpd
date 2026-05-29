# The full Hamiltonian $\hat H_{\rm full}$ — term by term

**Source:** Fig. 1c ("complete quantum description") of
[Minev et al., arXiv:2010.00620](https://arxiv.org/abs/2010.00620).

Cleaned up (with $\hbar=1$, so angular frequencies *are* energies, matching the figure):

$$
\hat H_{\mathrm{full}}
= \underbrace{\sum_{m}\left(\omega_m' - i\kappa_m\right)\hat a_m^\dagger \hat a_m}_{\text{linear: dressed + lossy modes}}
\;+\;
\underbrace{\sum_{\alpha,\beta} C_{\alpha,\beta}\prod_{m}\hat a_m^{\dagger\,\beta_m}\,\hat a_m^{\alpha_m}}_{\text{nonlinear: all orders}}
$$

Conceptually $\hat H_{\rm full}=\hat H_{\rm lin}+\hat H_{\rm nl}$: the first term is the
linearized Josephson circuit (plus loss), the second collects everything nonlinear.

---

## Term 1 — dressed, dissipative oscillators

$$\sum_{m}\left(\omega_m' - i\kappa_m\right)\hat a_m^\dagger \hat a_m$$

| Symbol | Meaning |
|---|---|
| $\sum_m$ | over the $M$ **eigenmodes of the *linearized* circuit** (each junction $\to$ small-signal inductor $L_j=\varphi_0^2/E_j$). These are the **hybridized** modes — qubit-like and cavity-like already mixed, so the linear coupling $g$ is *not* a separate term, it is absorbed here. |
| $\hat a_m^\dagger\hat a_m \equiv \hat n_m$ | photon-number operator of mode $m$. The term is just $\sum_m\omega_m'\hat n_m$: independent harmonic oscillators. |
| $\omega_m'$ | the **dressed** (observed) frequency, **not** the bare eigenfrequency $\omega_m$. Quantum fluctuations acting through the nonlinearity shift it by the **Lamb shift**: $\omega_m'\approx\omega_m-\Delta_m$, $\ \Delta_m=\tfrac12\sum_n\chi_{mn}$. $\omega_m$ is only an *intermediate* quantity. |
| $-i\kappa_m$ | **dissipation**, folded into the same term as a non-Hermitian part. $\kappa_m$ = mode energy-decay rate (linewidth) $\Rightarrow$ amplitude $\propto e^{-\kappa_m t}$. |

Loss is computed on equal footing from the **loss EPRs**:

$$
\frac{1}{Q_m}=\frac{\kappa_m}{\omega_m}=\sum_\ell \frac{p_{m\ell}}{Q_\ell},
$$

summing dissipative elements $\ell$ (bulk/surface dielectrics, seams, metals), each weighted
by its participation $p_{m\ell}$ and intrinsic quality $Q_\ell$.

---

## Term 2 — nonlinear interaction, to all orders

$$\sum_{\alpha,\beta} C_{\alpha,\beta}\prod_{m}\hat a_m^{\dagger\,\beta_m}\,\hat a_m^{\alpha_m}$$

A single compact term holding **every** nonlinear process the junctions mediate.

| Symbol | Meaning |
|---|---|
| $\alpha=(\alpha_1,\dots,\alpha_M)$ | **multi-index**: $\alpha_m$ = number of times mode $m$ is *annihilated* in a monomial. |
| $\beta=(\beta_1,\dots,\beta_M)$ | **multi-index**: $\beta_m$ = number of times mode $m$ is *created*. |
| $\prod_m \hat a_m^{\dagger\,\beta_m}\hat a_m^{\alpha_m}$ | a **normal-ordered** monomial (all $\hat a^\dagger$ to the left) — required for the correct treatment of vacuum fluctuations. |
| $C_{\alpha,\beta}$ | c-number coupling strengths, **fixed entirely by the EPRs**. |

The couplings follow from expanding each junction's nonlinear energy with the flux operator
written in the mode basis:

$$
\hat\varphi_j=\sum_m \varphi_{mj}\left(\hat a_m+\hat a_m^\dagger\right),
\qquad
\varphi_{mj}=s_{mj}\sqrt{\frac{p_{mj}\,\hbar\omega_m}{2E_j}} ,
$$

with EPR sign $s_{mj}=\pm1$. For a tunnel junction
$E_j^{\rm nl}(\hat\varphi_j)=-E_j\big[\cos\hat\varphi_j-1+\tfrac12\hat\varphi_j^2\big]$, so every
$C_{\alpha,\beta}$ is a sum of products of the $\varphi_{mj}$ — i.e. **products of $\sqrt{p_{mj}}$**.

**Constraints / selection rules**

- Hermiticity: $C_{\alpha,\beta}=C_{\beta,\alpha}^{*}$.
- For a pure cosine, only monomials with **even** total order $\sum_m(\alpha_m+\beta_m)$ survive (orders $4,6,\dots$; odd orders appear with flux bias / asymmetry).

---

## Reduction to the familiar Kerr Hamiltonian

Dispersive regime $\to$ keep only **excitation-number-conserving** monomials ($\alpha=\beta$) at
quartic order. Term 2 becomes Eq. (25):

$$
\hat H_{\rm eff}=\sum_m\big(\omega_m-\Delta_m\big)\hat a_m^\dagger\hat a_m
-\sum_m \frac{\alpha_m}{2}\,\hat a_m^{\dagger 2}\hat a_m^{2}
-\sum_{m<n}\chi_{mn}\,\hat a_m^\dagger\hat a_m\,\hat a_n^\dagger\hat a_n .
$$

For a single junction $j$ (Eq. 26):

$$
\chi_{mn}=\frac{\hbar\,\omega_m\omega_n}{4E_j}\,p_{mj}\,p_{nj},
\qquad
\alpha_m=\frac{\chi_{mm}}{2},
\qquad
\Delta_m=\tfrac12\sum_n\chi_{mn}.
$$

So $C_{\alpha=\beta}=-\alpha_m/2$ for the $\hat a_m^{\dagger2}\hat a_m^2$ monomial and $-\chi_{mn}$
for the $\hat n_m\hat n_n$ monomial.

---

## Where this lives in code (`qiskit-metal` / `pyEPR`)

- **Term 1 frequencies & losses** — eigenmode solve + `calc_energy_electric/magnetic`
  (`renderer_ansys/ansys_renderer.py`), loss EPRs from the `dissipatives` dict.
- **Term 2, all orders** — `pyEPR.QuantumAnalysis.analyze_all_variations(cos_trunc, fock_trunc)`
  keeps the full cosine numerically rather than truncating to the quartic $\chi$ formula.

**Punchline:** Term 1 says *where the energy sits and how fast it leaks*; Term 2 says *how the
modes talk* — and **all** of it ($\omega_m',\kappa_m,C_{\alpha,\beta}$) is set by the EPRs
$p_{mj}$, $p_{m\ell}$ and signs $s_{mj}$ from a single eigenmode solve.
