# Physical origin of the Lamb shift (one resonator + one qubit)

Self-contained. The name is borrowed from atomic physics; the analogy is exact, with a twist
specific to the (an)harmonic-oscillator setting.

> Rendering: display equations are images (Warp); inline math is compilable `$...$`.

## Notation

| Symbol | Meaning |
|---|---|
| $q,c$ | qubit / resonator (cavity) mode |
| $\hat a_m,\hat a_m^\dagger,\hat n_m$ | ladder / number operators of mode $m$ |
| $g$ | linear qubit–resonator coupling |
| $\omega_m^{\rm bare}$ | bare, $g=0$ uncoupled frequency |
| $\omega_m$ | linear eigenmode (normal-mode) frequency — includes $g$, excludes Lamb shift |
| $\Delta_m$ | Lamb shift of mode $m$ |
| $\omega_m^{\rm meas}=\omega_m-\Delta_m$ | dressed / measured frequency |
| $\alpha_m$ | anharmonicity = self-Kerr ($=\tfrac12\chi_{mm}$) |
| $\chi_{qc}$ | cross-Kerr / dispersive shift |

## 1. The original Lamb shift (atoms, 1947)

H $2S_{1/2}$ and $2P_{1/2}$, degenerate in Dirac theory, are split. Cause: the electron couples to
the **quantized EM vacuum** — it emits and reabsorbs *virtual* photons even with no real photons
present, jiggling it and shifting its levels. **Lamb shift = level shift from coupling to vacuum
fluctuations.**

## 2. The circuit-QED analog

The qubit frequency is shifted by the **zero-point fluctuations** of the fields it couples to — even
with **zero real photons** in the resonator, whose vacuum still carries a half-quantum:

![vacuum half-quantum](img/15_vacuum.svg)

Those fluctuations, felt through the nonlinear (Kerr) coupling, pull the qubit frequency — that is
$\Delta_q$.

## 3. The twist: it is a *nonlinear* effect

**Two *linearly* coupled *harmonic* oscillators have no Lamb shift** — bilinear coupling just makes
normal modes, whose exact frequencies are the EPR eigenmodes $\omega_m$. The $g$-renormalization is
*already inside* $\omega_m$. The Lamb shift appears only on restoring the **anharmonicity**:

![frequency ladder: linear vs nonlinear](img/15_ladder.svg)

No nonlinearity ⇒ no Kerr ⇒ $\Delta_m=0$. (For a *two-level* qubit the JC shift $g^2/\Delta$ is
called a Lamb shift, but only because a two-level system is intrinsically nonlinear; the oscillator
picture separates the linear piece cleanly into $\omega_m$.) This is the "$+g$" vs "$+$Lamb" arrows
of [doc 12](12-bare-vs-dressed-frequencies.md).

## 4. Where the number comes from — vacuum half-quanta

The Lamb shift is the Kerr/cross-Kerr acting on the $\tfrac12$-quantum of vacuum in every mode:

![Lamb shift decomposition](img/15_decomp.svg)

- **$\alpha_q=\tfrac12\chi_{qq}$** — qubit dressed by its **own** zero-point motion sampling its
  anharmonic potential (self-dressing).
- **$\tfrac12\chi_{qc}$** — qubit dressed by the **resonator's** vacuum via the cross-Kerr. *This is
  the direct analog of the atomic Lamb shift* (dressed by another mode's photonic vacuum).

The resonator likewise gets $\Delta_c=\tfrac12\chi_{cc}+\tfrac12\chi_{qc}$ (dominated by the
qubit-vacuum term, since $\alpha_c$ is tiny).

**Virtual-photon view:** equivalently, the qubit virtually creates/reabsorbs excitations through the
nonlinear terms; second-order energy denominators give $\Delta$. Same physics, PT language.

## Read in the paper

"$\Delta$… due to the dressing of this nonlinear mode by quantum fluctuations of the fields" (after
Eq. 8); $\Delta_m=\tfrac12\sum_n\chi_{mn}$ (below Eq. 25). See also [doc 05](05-chi-as-dispersive-shift.md)
(Lamb vs dispersive; the $2n+1$ split) and [doc 12](12-bare-vs-dressed-frequencies.md) (frequency ladder).
