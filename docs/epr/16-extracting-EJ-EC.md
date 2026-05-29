# Extracting $E_J$ and $E_C$ from an EPR analysis

Self-contained. The honest answer: EPR gives you only **one** genuinely bare parameter ($E_J$, the input); a clean $E_C$ comes from **electrostatics** (the capacitance matrix), not from the hybridized eigenmode spectrum.

> Rendering: display equations are images (Warp); inline math is compilable `$...$`.

## Notation

| Symbol | Meaning |
|---|---|
| $E_J$ | Josephson energy |
| $E_C=e^2/2C_\Sigma$ | charging energy; $C_\Sigma$ = total capacitance shunting the junction |
| $L_J=\varphi_0^2/E_J$ | junction linear inductance — the EPR simulation **input** |
| $\alpha_q$ | qubit anharmonicity (EPR `chi_ND[q,q]`) |
| $p_{qJ}$ | participation of the junction in the qubit mode |
| $\omega_q$ | dressed qubit mode frequency |
| $f_{\rm ND}$ | dressed qubit $0\to1$ frequency (EPR `f_ND[q]`) |

## The core tension

$(E_J,E_C)$ describe the **decoupled** qubit subcircuit (junction + its shunt capacitor). But EPR returns **hybridized eigenmodes** of the *coupled* system. So **every eigenmode quantity — frequencies *and* anharmonicities — carries the hybridization.** Only the input is bare.

In particular the anharmonicity is *not* a clean bare-qubit property:

![alpha_q is hybridized](img/16_alpha_hyb.svg)

It depends on the participation $p_{qJ}<1$ (the mode is partly cavity-like, diluting the anharmonicity) and on the dressed $\omega_q$. Deducing $E_C$ from $\alpha_q$ therefore folds in the hybridization.

## Where a clean $E_C$ comes from: electrostatics

$E_C$ is a capacitance, read from the **Maxwell capacitance matrix** — a circuit-element quantity that exists *before* any normal mode forms, independent of eigenmode hybridization:

![E_C from capacitance](img/16_ec_cap.svg)

In the qiskit-metal stack this is the **Q3D capacitance / LOM** path (`analyses/quantization/lom_core_analysis.py`) — a *separate electrostatic* analysis from EPR.

The two analyses are complementary:
- **EPR (eigenmode):** $E_J$ (input) + the hybridized, dressed Hamiltonian ($\omega,\alpha,\chi$).
- **Q3D / LOM (electrostatic):** $E_C=e^2/2C_\Sigma$, bare.

Clean combination: $E_J$ from the EPR input, $E_C$ from the capacitance matrix.

## What's bare vs tainted

| Quantity | Clean (bare) source | Tainted source |
|---|---|---|
| $E_J$ | EPR input $\varphi_0^2/L_J$ | — |
| $E_C$ | $e^2/2C_\Sigma$ from Q3D/LOM capacitance | $\alpha_q$ from EPR (hybridized, $\sim$few-% off) |

## If EPR is all you have

$\alpha_q$ is your best handle on $E_C$, but it is the hybridized-mode value: $\alpha_q\approx E_C\times(\text{participation dilution})$. For a transmon the qubit mode is junction-dominated ($p_{qJ}\approx0.98$–$0.999$), so the error is a few percent — usually acceptable, but it is a hybridization-tainted $E_C$, not a bare one. For precision, de-hybridize using $p_{qJ}$, or use the capacitance route.

If you do this, solve the **bare** CPB $\hat H=4E_C(\hat n-n_g)^2-E_J\cos\hat\phi$ (e.g. with `QPD(e_j_hz,e_c_hz)`) at the known $E_J$ for the $E_C$ that reproduces $\alpha_q$ — but read the result as "the $E_C$ consistent with the *measured/hybridized* anharmonicity," not the bare charging energy.

## Why not use the dressed frequency `f_ND`?

It is the *most* contaminated quantity — hybridized **and** Lamb-shifted. Note, though, that the bare CPB $f_{01}=\sqrt{8E_JE_C}-E_C$ already contains the dominant ($\sim$200 MHz) self-anharmonicity Lamb shift (the $-E_C$); the residual difference from `f_ND` is only the small cross-Kerr Lamb shift $\tfrac12\chi_{qc}$ plus the linear hybridization. So treat `f_ND` as a **consistency check**, never as the bare $f_{01}$.

## Note vs. the EPR paper

The paper (arXiv:2010.00620) takes $E_J$ as input and reports $\alpha,\chi$; it does not parametrize in $(E_J,E_C)$. The bare $E_C$ lives in the electrostatic/capacitance description (LOM), which is complementary to EPR and underlies the `qpd` `fit_quantum_capacitance` workflow.
