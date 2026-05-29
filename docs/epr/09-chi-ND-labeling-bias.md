# The labeling issue in χ_ND, and why it biases dispersive-shift estimates

See [07](07-chi-ND-vs-perturbation.md)/[08](08-H_ND-in-kerr-form.md) for what χ_ND is. Here: how the
**eigenstate labeling** step works, a worked example where it fails, and why the failure is a
*systematic* bias.

## How labeling works

Diagonalizing H_ND gives *dressed* eigenstates |ψ_k⟩. The finite-difference formulas need to know
which |ψ_k⟩ "is" |1,1⟩, |1,0⟩, etc. Standard rule (pyEPR): assign bare label (n_a, n_b) to the
eigenstate of **largest overlap** |⟨n_a,n_b|ψ_k⟩|².

- **Weak nonlinearity:** each dressed state ≈ 99% one bare state → unambiguous → χ_ND ≈ χ_O1. Fine.
- **Strong hybridization:** overlaps drop toward 0.5 → assignment becomes ambiguous → trouble.

## Worked example

Two modes: true cross-Kerr χ0 = 2 MHz, plus a cubic term g(â†²b̂ + â²b̂†) — a non-number-conserving
piece that lives in the full cosine but is *dropped* in the Kerr/`chi_O1` reduction. It couples
|1,1⟩ ↔ |3,0⟩, resonant when ω_b ≈ 2ω_a − 3α_a (here ≈ 9402 MHz). Sweeping ω_b:

![chi_ND bias near resonance](img/09_labeling.png)

| ω_b (MHz) | χ_ND (MHz) | overlap of \|1,1⟩ | eigenstate k(1,1) | k(3,0) |
|---|---|---|---|---|
| 9300 | 5.4 | 0.967 | 4 | 5 |
| 9395 | 18.1 | 0.588 | 4 | 5 |
| 9402 | 21.3 | 0.500 | 4 | 5 |
| 9410 | −14.3 | 0.600 | **5** | **4** |
| 9500 | −2.2 | 0.964 | 5 | 4 |

True value is **2 MHz throughout**, yet χ_ND swings from +21 to −14, and the eigenstate index for
|1,1⟩ **swaps** (4 ↔ 5) across the crossing.

## The 2×2 model

Bare |1,1⟩ vs spectator |s⟩ = |3,0⟩, detuning ε, coupling V = ⟨s|H_nl|1,1⟩ = g√6:

![avoided crossing energies](img/09_2x2.svg)

Only E_11 sits near the resonance (E_00, E_10, E_01 do not), so the finite difference inherits the
repulsion directly:

![chi_ND bias formula](img/09_bias.svg)

## Why it's a *systematic* bias (not noise)

1. **Sign-definite.** Repulsion always pushes labeled |1,1⟩ *away* from the spectator: −V²/ε has a
   fixed sign on each side of the resonance. Sweeping E_J/geometry/photon number varies it smoothly
   but one-sidedly — it does not average out.
2. **Contaminates only E_11.** Just one of the four levels is near resonance, so the shift enters
   χ_ND undiluted.
3. **Label swap at the crossing.** At ε = 0 the eigenstates are 50/50; max-overlap is a tie and
   picks a branch arbitrarily. For ε<0 vs ε>0 the "mostly-|1,1⟩" eigenstate is a *different index*
   (k: 4 ↔ 5 above). A mislabel grabs the wrong level → discontinuous jump / sign flip in χ_ND.

## Why it matters physically

This is not a dismissable artifact — it *is* the real breakdown of the dispersive picture
(measurement-induced state transitions; χ going photon-number-dependent). The genuine cross-Kerr is
conflated with an accidental resonance against a **non-readout** state. χ_ND returns a correct level
difference, but a **biased estimate of the dispersive shift**.

Contrast:
- `chi_O1` stays flat at 2 MHz (it dropped the cubic term) → **misses** the real repulsion.
- `chi_ND` → **over-attributes** the repulsion to dispersion.

Near such resonances neither equals "the dispersive shift"; you must identify the resonance.

**Practical tell:** if the max overlap of the labeled state isn't ≈ 1, distrust χ_ND. Check
`fock_trunc` convergence and look for a nearby bare state it's colliding with.
