# Extracting $(E_C, g, \omega_r)$ for a Transmon–Resonator System

A pedagogical note on the two simulation-based routes to the bare Hamiltonian
parameters, the physics behind each, the relevant qiskit-metal functions, and
how to interpret the (generally different) numbers they produce.

---

## 1. The goal and why it is not trivial

We want the parameters of the **bare qubit + bare resonator + coupling** Hamiltonian

$$
\frac{\hat H}{\hbar} = \sum_i \omega_{q,i} |i\rangle\langle i|\otimes \mathbb{1}_r +
\omega_r \mathbb{1}_q\otimes \hat a^\dagger \hat a +
g (\hat n - n_g)_q \otimes (\hat a + \hat a^\dagger),
$$

where the transmon is a Cooper-pair box (CPB)

$$
\hat H_{\rm CPB} = 4E_C(\hat n - n_g)^2 - E_J\cos\hat\varphi,
\qquad [\hat\varphi,\hat n]=i .
$$

The model has **four parameters**: $E_J, E_C, \omega_r, g$. We assume $E_J$ is known
(from the junction's normal-state resistance via Ambegaokar–Baratoff, or design
intent), and we want the remaining three: $(E_C, g, \omega_r)$.

**What we need.** The two dressed normal-mode frequencies $\Omega_\pm$ are the most
obvious observables, but they are not sufficient on their own to pin down all three
parameters — we also need an **independent third observable** that resolves the
qubit's charging scale. There are two natural places to get it: the **mode shapes**
(junction participations, from an eigenmode solve) or the **capacitance network**
(from an electrostatic solve). These define the two methods below.

### Linearized building blocks (used by both methods)

Expanding the transmon to harmonic order ($E_J \gg E_C$, phase localized near 0):

$$
\omega_q = \frac{\sqrt{8 E_C E_J}}{\hbar},\qquad
\varphi_{\rm zpf} = \Big(\frac{2E_C}{E_J}\Big)^{1/4},\qquad
n_{\rm zpf} = \frac{1}{2}\Big(\frac{E_J}{2E_C}\Big)^{1/4},
$$

with $\hat\varphi_J = \varphi_{\rm zpf}(\hat b + \hat b^\dagger)$ and
$\hat n = i n_{\rm zpf}(\hat b^\dagger - \hat b)$. Charge is the *momentum*
quadrature, which is why capacitive coupling lands on the antisymmetric
combination $\hat b^\dagger - \hat b$ rather than $\hat b + \hat b^\dagger$.
Here $\hat\varphi_J$ is the junction phase $\hat\varphi$ of the CPB, now written in
the harmonic representation (it picks up the higher powers of the cosine only when
the nonlinearity is restored).
The linearized coupling is $\tilde g (\hat b^\dagger-\hat b)(\hat a+\hat a^\dagger)$
with $\tilde g = g n_{\rm zpf}$, and the normal-mode frequencies are

$$
\Omega_\pm^2 = \frac{\omega_r^2+\omega_q^2}{2}
\pm \frac{1}{2}\sqrt{(\omega_r^2-\omega_q^2)^2 + 16 \tilde g^2 \omega_r\omega_q}
\qquad\text{(exact),}
$$

or, in the dispersive/RWA limit, the simpler hopping form

$$
\Omega_\pm = \frac{\omega_r+\omega_q}{2}
\pm \frac{1}{2}\sqrt{(\omega_r-\omega_q)^2 + 4\tilde g^2}.
$$

---

## 2. Method A — EPR (eigenmode) approach

**Data source:** a full-wave *eigenmode* simulation (Ansys HFSS) of the coupled
structure, post-processed by Energy-Participation-Ratio analysis. It returns, for
each normal mode $m$:

- the dressed mode frequency $\Omega_m$, and
- the **junction energy participation** $p_m$ = fraction of mode $m$'s inductive
  energy stored in the Josephson junction.

### The physics: the junction zero-point phase is *partitioned* among the modes

The junction lives on the bare qubit, so its phase distributes over the normal
modes as $\hat\varphi_J = \sum_m \varphi_{J,m}(\hat c_m + \hat c_m^\dagger)$. Two
relations connect this to what EPR reports:

**(i) Participation ↔ junction phase zpf (exact definition).** With junction
inductive energy $\tfrac{E_J}{2}\hat\varphi_J^2$ and mode inductive zpf energy
$\tfrac14\hbar\omega_m$,

$$
p_m = \frac{2 E_J \varphi_{J,m}^2}{\hbar\omega_m}
\quad\Longleftrightarrow\quad
\varphi_{J,m}^2 = \frac{p_m \hbar\omega_m}{2E_J}.
$$

**(ii) Conservation of the junction zpf (orthogonal/RWA limit).** The junction sits
entirely on the bare qubit, so the junction phase *is* the bare-qubit phase,
$\hat\varphi_J = \varphi_{\rm zpf}(\hat b + \hat b^\dagger)$. Its scale is the
transmon's zero-point fluctuation: from the harmonic CPB — an oscillator with
"mass" $1/8E_C$ and frequency $\omega_q = \sqrt{8E_CE_J}$ — the ground-state position
variance is

$$
\varphi_{\rm zpf}^2 = \langle 0|\hat\varphi^2|0\rangle = \frac{1}{2 m \omega_q} = \sqrt{\frac{2E_C}{E_J}}.
$$

Hybridization is a number-conserving (RWA) beam-splitter, so the bare→normal map is
an **orthogonal** rotation $\hat b = \sum_m u_m\hat c_m$ with $\sum_m u_m^2 = 1$. It
merely *redistributes* the bare-qubit phase among the dressed modes,
$\varphi_{J,m} = \varphi_{\rm zpf} u_m$, and an orthogonal rotation preserves the
total — so the junction zpf is conserved:

$$
\sum_m \varphi_{J,m}^2 = \varphi_{\rm zpf}^2 \sum_m u_m^2 = \varphi_{\rm zpf}^2 = \sqrt{\frac{2E_C}{E_J}}.
$$

With counter-rotating terms the map is symplectic and the ground state is squeezed,
so this conservation acquires an $O((\tilde g/\omega)^2)$ correction; the exact,
RWA-free replacement is derived in Appendix A.

### The third equation: the participation sum rule

Combining (i) and (ii) gives a remarkably clean result — **the bare qubit frequency
is the participation-weighted average of the dressed frequencies:**

$$
\boxed{\omega_q = \sum_m p_m \Omega_m,
\qquad
E_C = \frac{(\hbar\omega_q)^2}{8 E_J} = \frac{\hbar^2\big(\sum_m p_m\Omega_m\big)^2}{8 E_J}}
$$

This single equation is what the linear spectrum lacked. The $p_m$ supply the
absolute junction zpf scale — i.e. $E_C$ — that the eigenfrequencies alone cannot.

### Finishing $\omega_r$ and $g$

With $\omega_q$ known, the remaining two parameters drop out of $\Omega_\pm$:

$$
\omega_r = \Omega_+ + \Omega_- - \omega_q
\quad\text{(RWA trace; exact: } \omega_r^2 = \Omega_+^2+\Omega_-^2-\omega_q^2),
$$

$$
\tilde g = \tfrac12\sqrt{(\Omega_+-\Omega_-)^2 - (\omega_q-\omega_r)^2},
\qquad g = \frac{\tilde g}{n_{\rm zpf}}.
$$

### Built-in consistency check

You measured four numbers $(\Omega_\pm, p_\pm)$ but used three combinations. The
fourth checks the 2-mode model: the participation **ratio** must agree with the
mixing angle implied by $(\omega_q,\omega_r,\tilde g)$,

$$
\frac{p_{\rm qubit\text{-}like} \Omega_{\rm qubit\text{-}like}}{p_{\rm res\text{-}like} \Omega_{\rm res\text{-}like}} = \cot^2\theta\Big|_{\tan 2\theta = 2\tilde g/(\omega_q-\omega_r)}.
$$

Disagreement signals that higher modes participate (the 2-mode reduction is leaking).

### Equivalent variants (same data, same answer)

- **$\phi^4$ coefficient matching** — expand both the lumped model and the
  HFSS+nonlinearity Hamiltonian to quartic order and match coefficients. Equivalent
  to the participation route for fixing the linear parameters.
- **Full `cos_trunc` numerical diagonalization** — keep the cosine to high order and
  diagonalize in the Fock basis; this is what you use for the *anharmonic* physics
  (anharmonicity, cross-Kerr $\chi$) once the linear parameters are set.

### qiskit-metal

```python
from qiskit_metal.analyses.quantization import EPRanalysis   # energy_participation_ratio.py
# Eigenmode simulation backend:
#   qiskit_metal.analyses.simulation.EigenmodeSim  (Ansys HFSS eigenmode)
# Spectrum/anharmonicity is delegated to pyEPR:
#   renderer.epr_spectrum_analysis(cos_trunc=8, fock_trunc=7)
# Key knobs:
#   cos_trunc  -> highest Taylor order of cos(phi): keeps phi^4 ... phi^(2*cos_trunc)
#   fock_trunc -> Fock-space size per mode
```

Relevant class: `EPRanalysis` (`qiskit_metal/analyses/quantization/energy_participation_ratio.py`),
backed by `EigenmodeSim` (`.../simulation/eigenmode.py`).

---

## 3. Method B — LOM (capacitance) approach

**Data source:** an *electrostatic* capacitance-matrix simulation (Ansys Q3D, or
the open-source Elmer/Gmsh path) plus the junction inductance $L_J$ and a
transmission-line model for the resonator.

### The physics: electrostatics resolves what the spectrum hides

The Maxwell capacitance matrix separates the **self-capacitance** of the qubit node
from the **coupling capacitance** — two distinct entries that the linear spectrum
only ever sees in combination. That separation is the missing third equation, but
sourced statically rather than dynamically.

$$
E_C = \frac{e^2}{2 C_\Sigma}
\qquad\text{(qubit-node self-capacitance)} .
$$

The capacitive coupling follows from the off-diagonal (coupling) capacitance and the
modes' zero-point voltages (Koch et al. form):

$$
\hbar g = 2 e \beta V^0_{\rm rms} n_{\rm zpf},
\qquad
\beta = \frac{C_c}{C_\Sigma},
\qquad
V^0_{\rm rms} = \sqrt{\frac{\hbar\omega_r}{2 C_r}} .
$$

In matrix language this is $g_{ij}\propto [C^{-1}]_{ij} Q_{{\rm zpf},i}Q_{{\rm zpf},j}$ —
the off-diagonal of the **inverse** capacitance matrix weighted by zero-point charges.

### Where $\omega_r$ comes from (and why it is *bare*, not dressed)

LOM does **not** diagonalize the coupled system to get $\omega_r$. The capacitance
matrix is purely electrostatic — the junction inductance $L_J$ never enters it — so
the resonator frequency is built from the resonator's own $L_r$ and its diagonal
self-capacitance (loaded by $C_c$ but **not** hybridized with the qubit mode):

$$
\omega_r = \frac{1}{\sqrt{L_r (C_r + C_c)}} .
$$

This captures *capacitive loading* (effect of the qubit metal) without *mode
hybridization* (which requires $L_J$). It is therefore the bare $\omega_r$ that
belongs in the model Hamiltonian, with the hybridization held separately in $g$.

### qiskit-metal

```python
from qiskit_metal.analyses.quantization import LOManalysis    # lumped_oscillator_model.py
# Capacitance-extraction backend:
#   qiskit_metal.analyses.simulation.LumpedElementsSim   (Ansys Q3D / Elmer)
# Core conversions (from pyEPR.calcs.convert.Convert):
#   EC = Convert.Ec_from_Cs(1 / c_inv_k[ss, ss])     # E_C from self-capacitance
#   EJ = Convert.Ej_from_Lj(1 / l_inv_k[ss, ss])     # E_J from junction inductance
#   gs = analysis.compute_gs(CouplingType.CAPACITIVE)  # couplings from C-matrix
# Also returns the cross-Kerr matrix:  ham_res['chi_in_MHz']
```

Relevant class: `LOManalysis` (`qiskit_metal/analyses/quantization/lumped_oscillator_model.py`),
core logic in `lom_core_analysis.py`, backed by `LumpedElementsSim`
(`.../simulation/lumped_elements.py`).

---

## 4. Interpreting different values from the two methods

**In the ideal limit they are provably equal.** If the device were *exactly* a
two-mode lumped circuit (one resonator mode, one junction-bearing transmon mode,
bilinear capacitive coupling), both methods recover the *same* $(E_C,g,\omega_r)$.
This is forced by the transmon having a single $E_C$ that simultaneously sets the
capacitance (LOM) and the junction zpf / anharmonicity (EPR).

**Which is more fundamental? — EPR.** Run to full convergence, EPR is the more
fundamental linear extraction. It solves the true full-wave linear Maxwell problem
and reads off the exact spectrum $\Omega_\pm$ and participations $p_\pm$; all of its
approximations — the number of modes kept in $\omega_q = \sum_m p_m\Omega_m$, the
mesh, and the (removable) RWA in the sum rule — are *convergence-controllable* and
vanish with effort. LOM instead builds a lumped LC + transmission-line circuit from a
*quasistatic* capacitance solve, and that quasistatic/lumped reduction is a
**structural** approximation that does **not** vanish with mesh refinement. Note that
the 2-mode reduction is shared by both methods, and the cosine truncation is not used
in the linear participation route — so neither makes EPR rougher. In principle EPR is
at least as faithful as LOM, and more so.

**Then why ever use LOM? — purely practical.** EPR's fidelity must be paid for: a
full-wave eigenmode solve is expensive, and you need *many* modes, because the
high-frequency tail of $\omega_q = \sum_m p_m\Omega_m$ converges slowly — and that
tail is exactly what fixes $E_C$. An under-converged EPR run therefore *underestimates*
$E_C$. LOM's electrostatic $E_C = e^2/2C_\Sigma$ captures all electrostatic
contributions in one cheap, well-converged solve. So the reason to reach for LOM is
**cost and convergence — most acutely for $E_C$ — not any fundamental superiority.**

**Convention — how to read a disagreement, not a tiebreaker.** The bare parameters
$\omega_r$ and $g$ are not physical observables; only the spectrum and participations
are. Both methods fix a bare/dressed convention to repackage the same physics, and
those conventions need not be identical, so part of any LOM-vs-EPR gap is convention
rather than error. This affects *interpretation*; because both methods share it, it
favors neither.

**How to act on a discrepancy:**

- *Small*: your 2-mode lumped Hamiltonian is faithful and EPR is well converged. Use
  either set.
- *Large*: either EPR has not converged (add modes, refine the mesh — check the $E_C$
  tail first) or the device has structure beyond the 2-mode model (extra
  package/slotline modes, junction non-ideality). Rule out EPR convergence before
  reading the residual as genuine model incompleteness.

> **Rule of thumb:** run *both* and compare. Agreement bounds *both* your model error
> *and* your EPR convergence; disagreement tells you to add modes or enrich the model.
> The cross-check is the value — it is not an argument for preferring one method.

---

## 5. Assumptions side-by-side

| Aspect | EPR (eigenmode) | LOM (capacitance) |
|---|---|---|
| **Primary simulation** | Full-wave eigenmode (HFSS) | Electrostatic capacitance matrix (Q3D / Elmer) |
| **Raw observables** | $\Omega_\pm$ and junction participations $p_\pm$ | Maxwell capacitance matrix $C_{ij}$ (+ $L_J$) |
| **Third equation for $E_C$** | Participation sum rule $\omega_q=\sum_m p_m\Omega_m$ | Self-capacitance $E_C=e^2/2C_\Sigma$ |
| **Nature of $E_C$** | Dynamic (AC, mode-frequency response) | Static (DC capacitance) |
| **Source of $\omega_r$** | Inversion of $\Omega_\pm$ (after $\omega_q$ fixed) | Resonator $L_r$ + diagonal $C$ (TL model) |
| **Source of $g$** | Mode gap $\tilde g=\tfrac12\sqrt{(\Omega_+-\Omega_-)^2-(\omega_q-\omega_r)^2}$ | Off-diagonal $C_c$ + zero-point voltages |
| **Linear coupling $g$** | Treated **exactly** (eigenmodes already hybridized) | Treated **exactly** (full $C$-matrix) |
| **Multimode physics** | Captured exactly, but the high-mode tail (sets $E_C$) converges slowly | Electrostatics all-mode; resonator via single TL mode |
| **Distributed fields** | Native (full-wave) | Approximated (lumped + TL) |
| **Key approximation** | Convergence-controllable: mode count in sum rule, mesh, removable RWA | Structural: quasistatic + lumped/TL reduction (does not vanish with mesh) |
| **In-principle fidelity** | Higher (full-wave; approximations vanish on convergence) | Lower (quasistatic/lumped error is irreducible) |
| **Anharmonicity / $\chi$** | Yes — via `cos_trunc` diagonalization | Yes — via Kerr from $C$-matrix (`chi_in_MHz`) |
| **Convention dependence** | Bare = junction-bearing mode (set by participations) | Bare = node basis of the $C$-matrix |
| **qiskit-metal class** | `EPRanalysis` (+ `EigenmodeSim`) | `LOManalysis` (+ `LumpedElementsSim`) |
| **Best trusted for** | Most faithful when converged; spectroscopic physics; distributed/multimode regimes | Cheap, robust, well-converged — especially $E_C$; fast design iteration |

---

## 6. Practical recipe

1. **LOM first** (cheap, electrostatic): get a robust $E_C$ and an initial
   $(\omega_r, g)$ for design iteration.
2. **EPR second** (full-wave): get $\Omega_\pm$ and $p_\pm$; extract
   $(E_C,\omega_r,g)$ via the participation sum rule, and the anharmonic physics via
   `cos_trunc` diagonalization.
3. **Compare** the two $(E_C,g,\omega_r)$ triples. Use the agreement as validation;
   read a discrepancy as a combined EPR-convergence and model-error bar (check the
   EPR mode count first).
4. **Use the chosen triple** in the bare-basis Hamiltonian of Section 1 (with the
   full cosine, not the linearized form) to compute whatever physics you need.

---

### Symbol reference

| Symbol | Meaning |
|---|---|
| $E_J, E_C$ | Josephson / charging energy |
| $\omega_q$ | Bare transmon (plasma) frequency, $\sqrt{8E_CE_J}/\hbar$ |
| $\omega_r$ | Bare resonator frequency |
| $g, \tilde g$ | Bare coupling; linearized coupling $\tilde g = g n_{\rm zpf}$ |
| $\Omega_\pm$ | Dressed (normal-mode) frequencies from the coupled system |
| $p_m$ | Junction energy participation of mode $m$ |
| $\varphi_{\rm zpf}, n_{\rm zpf}$ | Zero-point phase / charge fluctuations of the bare transmon |
| $C_\Sigma, C_c, C_r$ | Qubit self-, coupling, resonator capacitance |

---

## Appendix A. RWA-free (exact) extraction

The Section-2 recipe used the RWA in one place only: the sum rule
$\sum_m p_m\Omega_m=\omega_q$, which rested on the junction-zpf *conservation*
$\sum_m\varphi_{J,m}^2=\varphi_{\rm zpf}^2$ (valid for an orthogonal, i.e. RWA,
bare-to-normal map). The two Vieta relations were already exact. Here we replace
that one step by its exact counterpart and obtain a fully RWA-free, closed-form
extraction. (Every formula below was checked against exact Fock-space
diagonalization and an honest capacitively-coupled circuit, agreeing to machine
precision from weak coupling through strong coupling.)

### A.1 Setup

Set $\hbar=1$. The exact linearized Hamiltonian in quadratures is

$$
H=\frac{\omega_r}{2}(x_a^2+p_a^2)+\frac{\omega_q}{2}(x_b^2+p_b^2)+2\tilde g x_a p_b,
$$

with the junction phase on the bare qubit, $\hat\varphi_J=\varphi_{\rm zpf}(\hat b+\hat b^\dagger)=\sqrt2 \varphi_{\rm zpf} x_b$, so

$$
\langle\hat\varphi_J^2\rangle=2\varphi_{\rm zpf}^2\langle x_b^2\rangle .
$$

Two exact facts carry over from Section 2 and from the participation definition:

$$
\Omega_+^2+\Omega_-^2=\omega_r^2+\omega_q^2,\qquad
\Omega_+^2\Omega_-^2=\omega_r^2\omega_q^2-4\tilde g^2\omega_r\omega_q\quad\text{(Vieta)},
$$

$$
\varphi_{J,m}^2=\frac{p_m\Omega_m}{2E_J},\qquad
\sum_m\varphi_{J,m}^2=\langle\hat\varphi_J^2\rangle .
$$

The participation identity $\varphi_{J,m}^2=p_m\Omega_m/2E_J$ is exact per mode — it
follows from mode equipartition (total inductive energy $=\Omega_m/4$), which holds
for genuine normal modes and was confirmed numerically against the EPR energy-ratio
definition, **with no RWA assumed**. The only quantity that was RWA-approximated
before is the total junction-phase variance $\langle\hat\varphi_J^2\rangle$, which we
now compute exactly.

### A.2 Exact junction-phase variance

The normal-mode eigenvectors of the equations of motion (mode $\propto e^{-i\Omega_m t}$, normalized later) are, taking $p_b=1$,

$$
x_a=-\frac{2\tilde g \omega_r}{\omega_r^2-\Omega_m^2},\quad
p_a=\frac{2i\tilde g \Omega_m}{\omega_r^2-\Omega_m^2},\quad
x_b=\frac{i\Omega_m}{\omega_q},\quad p_b=1 .
$$

The symplectic norm $N_m=\mathbf v^\dagger(i \Omega_{\rm sympl})\mathbf v$ is

$$
N_m=2\Omega_m\left[\frac{4\tilde g^2\omega_r}{(\omega_r^2-\Omega_m^2)^2}+\frac{1}{\omega_q}\right],
$$

so each mode contributes $|x_b|^2/N_m$ to $\langle x_b^2\rangle$. Using the
determinant identity $(\omega_r^2-\Omega_m^2)(\omega_q^2-\Omega_m^2)=4\tilde g^2\omega_r\omega_q$ and the Vieta sum, the per-mode term collapses to

$$
\langle x_b^2\rangle=\sum_m\frac{|x_b|^2}{N_m},\qquad
\frac{|x_b|^2}{N_m}=\frac{\Omega_m(\omega_r^2-\Omega_m^2)}{2\omega_q(\Omega_+^2+\Omega_-^2-2\Omega_m^2)} .
$$

Summing the two modes gives the closed form

$$
\boxed{\langle x_b^2\rangle=\frac{\omega_q^2+\Omega_+\Omega_-}{2\omega_q(\Omega_++\Omega_-)}}
$$

which equals $\tfrac12$ in the decoupled limit, as it must.

### A.3 Exact third equation

Combining $\sum_m p_m\Omega_m=2E_J\langle\hat\varphi_J^2\rangle=2\omega_q\langle x_b^2\rangle$ (using $2E_J\varphi_{\rm zpf}^2=\omega_q$) with the boxed variance:

$$
\boxed{\sum_m p_m\Omega_m=\frac{\omega_q^2+\Omega_+\Omega_-}{\Omega_++\Omega_-}}
$$

This is the exact, RWA-free replacement for $\sum_m p_m\Omega_m=\omega_q$. In the
dispersive limit $\Omega_+\Omega_-\to\omega_q\omega_r$, $\Omega_++\Omega_-\to\omega_q+\omega_r$, the RHS returns $\omega_q$.

### A.4 Closed-form RWA-free extraction

Three exact equations in $(\omega_q,\omega_r,\tilde g)$, solved in sequence:

$$
\omega_q^2=(\Omega_++\Omega_-)(p_+\Omega_++p_-\Omega_-)-\Omega_+\Omega_-,\qquad E_C=\frac{\omega_q^2}{8E_J},
$$

$$
\omega_r^2=\Omega_+^2+\Omega_-^2-\omega_q^2,
$$

$$
\tilde g^2=\frac{\omega_r^2\omega_q^2-\Omega_+^2\Omega_-^2}{4 \omega_r\omega_q},\qquad g=\frac{\tilde g}{n_{\rm zpf}} .
$$

Every input is an observable ($\Omega_\pm$, $p_\pm$); no RWA enters. The recovered
bare parameters sit in the charge-basis convention — the split built from the
diagonal of $C^{-1}$ — which is the same one the Section-2 formulas converge to in
the dispersive limit.

### A.5 What is still approximate

Removing the RWA leaves two approximations that are inherent to the *target model*,
not to RWA: the **2-mode truncation** (extend by summing over all modes in
$\langle x_b^2\rangle$ and a many-mode Vieta) and the **linearization of the junction**
(the cosine's higher orders are the separate, nonlinear story of Section 2).
