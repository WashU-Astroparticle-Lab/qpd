# Writing H_ND with α_m and χ_mn (the Kerr form)

Starting Hamiltonian (full cosine; see [07](07-chi-ND-vs-perturbation.md)):

> H_ND = Σ_m ħω_m â_m†â_m − Σ_j E_j [cos φ̂_j − 1 + ½ φ̂_j²]

**You cannot write H_ND *exactly* with only α_m and χ_mn** — those capture only the leading quartic,
number-conserving part. That leading reduction is the Kerr Hamiltonian (= Eq. 25).

## Reduction

Expand the cosine and substitute φ̂_j = Σ_m φ_mj(â_m + â_m†):

![cosine expansion to quartic](img/08_expand.svg)

Normal-order and keep the number-conserving monomials (â_m†²â_m² and â_m†â_m â_n†â_n):

![Kerr-form Hamiltonian](img/08_kerrH.svg)

with

![alpha definition](img/08_alpha.svg)

![chi definition](img/08_chi.svg)

![Lamb shift](img/08_lamb.svg)

(Combinatorics: â†²â² enters φ̂⁴ with weight 6 → α_m; the mixed n̂_m n̂_n term with weight 6×4 = 24
→ χ_mn; the leftover â†â pieces give the Lamb shift Δ_m. Note α_m = ½ χ_mm.)

In EPRs: φ_mj² = p_mj ħω_m / (2E_j), which is what turns the φ-forms into the p-forms above
(matching Eq. 26 for χ and Eq. 9/10 for α).

## Why "≃", not "="

Writing H_ND with only α_m, χ_mn discards:

1. **higher-order Kerr** — the φ̂⁶ term → â†³â³, etc.;
2. **non-number-conserving / counter-rotating** terms — kept in H_ND, dropped here.

This α,χ form **is** the perturbative `chi_O1` Hamiltonian. `chi_ND` diagonalizes the full H_ND
*before* this reduction — exactly the `chi_O1` vs `chi_ND` distinction (see
[07](07-chi-ND-vs-perturbation.md)).

## Read in the paper

The boxed Hamiltonian is **Eq. (25)**; χ_mn is **Eq. (26)**; α_m is **Eq. (9/10)**; Δ_m = ½ Σ_n χ_mn
is stated below Eq. (25). It is the leading-order face of H_ND, i.e. of Eq. (17).

---

## Are α_m, χ_mn "redefined" between Eq. (9–11) and χ_ND?

No — there is **one definition**, evaluated on **two spectra**. The invariant, physical definition
is *spectral*: anharmonicity = curvature of the ladder; cross-Kerr = mixed second difference of the
energies E(n⃗). Here e_m = occupation vector with one quantum in mode m (m-th unit vector, a state
label):

![occupation-vector notation](img/legend_em.svg)

![one spectral definition](img/08_onedef.svg)

Feed it two different spectra:

![same definition, two spectra](img/08_twospectra.svg)

- **Perturbative spectrum** (eigenvalues of the Eq.-8 Kerr Hamiltonian) → recovers exactly the
  Eq. (9–11) formulas → `chi_O1`.
- **Exact spectrum** (numerical diagonalization of full H_ND) → `chi_ND`.

Same finite-difference operation, different input energies. So α_m, χ_mn are the **same observable**,
not a redefinition: Eqs. (9–11) are its leading-order analytic value, `chi_ND` its exact numerical
value. They agree in the weakly-nonlinear limit and diverge when higher-order terms make the true
spectrum non-Kerr.

**The subtlety you sensed:** because the exact spectrum is *not* perfectly Kerr, the finite-difference
α_m, χ_mn depend on *which* levels you pick (by convention the bottom, 0,1,2). The per-photon shift
at high n differs from χ_mn extracted at 0→1; capturing that needs *extra* parameters (6th-order
Kerr, …). That residual is exactly why H_ND can't be written with α_m, χ_mn alone — they are the
**leading spectral coefficients**, not exact Hamiltonian parameters.

> One definition, two evaluations. `chi_ND` is the faithful value of the same α_m, χ_mn; Eqs. (9–11)
> are its first-order approximation. pyEPR reports both because they are the same quantity computed
> two ways — and their gap is the breakdown gauge.

## Provenance: which definition is actually in the paper?

Important caveat on the two "definitions" above:

- **Perturbative / formula definition — IN the paper.** α_m, χ_mn are defined as the coefficients
  of the effective excitation-number-conserving Hamiltonian: **Eq. (25)** + the sentence right after
  it (general), and **Eq. (8)** (simple example). Their values are **Eq. (26)** (with α_m = χ_mm/2,
  Δ_m = ½ Σ_n χ_mn) and **Eqs. (9)–(12)** (example). Derivation: **Supplementary Section B**.

- **Spectral finite-difference definition (E_11 − E_10 − E_01 + E_00) — NOT in the paper.** The paper
  only states H_full "can be analytically or numerically diagonalized using various computational
  techniques [20]" and that pyEPR [95] extracts "its quantum spectrum," without writing the
  extraction formula. The explicit finite-difference recipe is **pyEPR's implementation convention**
  (and standard cQED practice from the black-box-quantization lineage, e.g. Nigg et al. 2012 =
  paper ref [4]).

So the "one definition, two spectra" framing is the *conceptually* unifying view, but only the
perturbative branch is written in the paper; the spectral branch lives in pyEPR.
