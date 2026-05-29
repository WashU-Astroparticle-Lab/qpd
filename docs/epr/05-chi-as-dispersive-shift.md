# How to see $\chi_{qc}$ is the dispersive shift, from Eq. (8)

Eq. (8):

> $\hat H_\text{eff} = (\omega_q - \Delta_q) \hat n_q + (\omega_c - \Delta_c) \hat n_c - \chi_{qc} \hat n_q \hat n_c$
>         $- \tfrac12 \alpha_q \hat n_q(\hat n_q - 1) - \tfrac12 \alpha_c \hat n_c(\hat n_c - 1)$

**Trick:** group the cross-Kerr term $-\chi_{qc} \hat n_q \hat n_c$ with the cavity number operator.

![dispersive grouping](img/05_dispersive.svg)

The cavity's effective frequency becomes $\omega_c - \Delta_c - \chi_{qc}\cdot\hat n_q$ — i.e. it is pulled by **$-\chi_{qc}$ per qubit excitation**:

- qubit in $|0\rangle$  $\to$  cavity at $\omega_c - \Delta_c$
- qubit in $|1\rangle$  $\to$  cavity at $\omega_c - \Delta_c - \chi_{qc}$

That qubit-state-dependent cavity frequency **is** the dispersive shift — the basis of dispersive readout. By symmetry of the $\hat n_q \hat n_c$ term, the qubit frequency likewise shifts by $-\chi_{qc}$ per cavity photon (the AC-Stark / number-splitting view).

## Read in the paper

The sentence introducing Eq. (8) already names $\chi_{qc}$ the "qubit–cavity dispersive shift (cross-Kerr coupling)." The grouping above is just how the $\hat n_q \hat n_c$ term realizes that label.

---

## Lamb shift vs. dispersive shift (don't confuse $\Delta_q$ and $\chi_{qc}$)

Both appear in Eq. (8). The clean distinction is **offset vs. slope** in the *partner* mode's photon number. The qubit $0\to1$ transition at fixed cavity occupation $n_c$ (from Eq. 8; the $\alpha_q$ term drops out for the $0\to1$ transition) is:

![Lamb = intercept, dispersive = slope](img/05_lamb_vs_disp.svg)

| | **Lamb shift $\Delta_q$** | **Dispersive shift $\chi_{qc}$** |
|---|---|---|
| What | vacuum renormalization of the frequency | state-dependent pull per partner excitation |
| In Eq. (8) | shifts $\omega_q \to \omega_q - \Delta_q$ | the $\hat n_q \hat n_c$ cross-Kerr term |
| Photon dependence | **none** (present even at $n_c = 0$) | **linear in $n_c$** (slope) |
| Origin | zero-point fluctuations dressing the nonlinearity | real excitations in the other mode |
| Observable | absolute dressed frequency | qubit-state $\leftrightarrow$ cavity-freq correlation (readout) |

They are two faces of the *same* nonlinearity. In fact the Lamb shift is built from the cross/self-Kerr acting on the half-quantum of vacuum in every mode:

![Lamb shift as sum of Kerr on vacuum](img/05_lamb_sum.svg)

i.e. each mode $n$ contributes $\tfrac12\chi_{mn}$ from its zero-point half-photon — that is *why* the vacuum ($n_c = 0$) frequency is already shifted.

**Read in the paper:** the text right after Eq. (8) ("$\Delta$… due to the dressing of this nonlinear mode by quantum fluctuations of the fields"); the experimental relation $\Delta_q = \alpha_q - \chi_{qc}/2$; and the general $\Delta_m = \tfrac12 \sum_n \chi_{mn}$ below Eq. (25)–(26).

---

## Aren't *both* shifts due to vacuum fluctuations?

Tempting, but **no — only the Lamb shift is a vacuum effect.** They are the *same* cross-Kerr interaction, evaluated on different things. Write it symmetrically (since $(\hat a+\hat a^\dagger)^2 \to 2\hat n+1$):

![symmetric cross-Kerr decomposition](img/05_symmetric.svg)

Equivalently, the qubit frequency is pulled by $\tfrac12\chi_{qc}\cdot(2n_c + 1)$:

![qubit freq (2n+1) decomposition](img/05_2n1.svg)

- the **"1"** is the cavity's zero-point half-quantum $\to$ **Lamb shift** (survives at $n_c = 0$);
- the **"$2n_c$"** is **real photons** $\to$ **dispersive / AC-Stark shift** (zero if the cavity is empty).

So:

- **Lamb shift** — yes, genuinely vacuum: the qubit dressed by the *zero-point* fluctuations of the cavity (and itself). It is the cross/self-Kerr acting on the $\tfrac12$'s, which is exactly why $\Delta_m = \tfrac12 \sum_n \chi_{mn}$.
- **Dispersive shift** — **not** vacuum. It is the response to *real* excitations $n_c$. Empty cavity $\Rightarrow$ no dispersive shift, only the Lamb shift.

**The grain of truth:** the *coefficient* $\chi_{qc}$ is itself set by the zero-point fluctuation *amplitudes* $\varphi_q$, $\varphi_c$ ($\chi_{qc} = \hbar\omega_q\omega_c p_q p_c / 4E_J$, and $\varphi_m^2 \propto$ ZPF). So vacuum sets the *strength* of the knob; real photons are what *turn* it. Strength (vacuum) vs. drive (real photons) — that's the clean split.
