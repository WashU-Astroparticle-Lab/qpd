# Do Eq. (11) and Eq. (12) contradict each other?

**No** — Eq. (12) is Eq. (11) re-expressed via Eqs. (9)–(10). They are the same number.

The relevant equations (Section II.A):

- Eq. (9):  $\alpha_q = \tfrac12 \chi_{qq} = p_q^2 \hbar\omega_q^2 / (8 E_J)$
- Eq. (10): $\alpha_c = \tfrac12 \chi_{cc} = p_c^2 \hbar\omega_c^2 / (8 E_J)$
- Eq. (11): $\chi_{qc} = p_q p_c \hbar \omega_q \omega_c / (4 E_J)$
- Eq. (12): $\chi_{qc} = \sqrt{\chi_{qq} \chi_{cc}} = 2\sqrt{\alpha_q \alpha_c}$

## One-line consistency check

Using $\chi_{mm} = 2\alpha_m$:

![consistency check](img/04_consistency.svg)

So (11) and (12) are identical.

## Why it *looks* like it might overconstrain

For a **single junction** shared by both modes, the EPR "vector" is rank-1, so the cross-Kerr **saturates** the Cauchy–Schwarz bound — equality, hence $\chi_{qc} = \sqrt{\chi_{qq} \chi_{cc}}$ *exactly*.

With **multiple junctions**, the cross-Kerr becomes a scalar product of EPR vectors (general **Eq. (26)–(28)**), and only the inequality $\chi_{mn} \le \sqrt{\chi_{mm} \chi_{nn}}$ holds.

## Read in the paper

- Derivation: **Supplementary Section B**.
- General matrix / scalar-product form: **Eqs. (26)–(28)**.
- The text between Eq. (11) and (13) explicitly notes $\chi_{qc}$ and $\alpha_q$ are *interdependent* (a single EPR $p_m$ fixes each mode's nonlinearity) — i.e. by design, not contradiction.
