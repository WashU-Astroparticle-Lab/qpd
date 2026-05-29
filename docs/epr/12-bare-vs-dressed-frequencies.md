# How to get the "bare" qubit and resonator frequencies (and why they aren't ω_q, ω_c)

"Bare" is overloaded. There are **three** frequencies, and the paper's ω_q, ω_c sit in the *middle*:

![frequency ladder](img/12_ladder.svg)

| Frequency | What it is | How to get it |
|---|---|---|
| **Linear eigenmode** ω_q, ω_c | junction → linear inductor L_J; **includes** g-hybridization, **excludes** Lamb shift | FE eigenmode solver directly (pyEPR `f_0`) |
| **Dressed / measured** ω_m − Δ_m | the physical transition you measure | subtract Lamb shift; pyEPR `f_1` (PT) or `f_ND` (numerical) |
| **Bare uncoupled** ω_m^bare | mode if g were off (isolated transmon, empty cavity) | **not** a native EPR output (see (c)) |

## (a) You want the linear-circuit mode frequencies

Those *are* ω_q, ω_c — read straight off the eigenmode simulation (`f_0`). But the paper's caveat:
they "will be significantly perturbed by the Lamb shifts Δ_m, and should be seen as an
**intermediate parameter**" (after Eq. 16/17). Not what you measure.

## (b) You want the physical, measurable frequency

Then *not* ω_q, ω_c. It is the dressed frequency:

![dressed = linear − Lamb shift](img/12_dressed.svg)

Get it by adding the Lamb shift (pyEPR `f_ND` / `f_1`), or just measure it. The paper compares to
experiment using exactly these "dressed mode frequencies ω_m − Δ_m" (Comparison section / Tables).

## (c) You want the truly bare, g-uncoupled frequencies

Also **not** ω_q, ω_c — and EPR does not hand these to you. The paper: the linear coupling g "is
fully factored in our analysis, and is implicitly handled in the extraction of the operators." The
eigensolver returns the *already-hybridized* modes. To recover g-uncoupled values, step outside the
standard flow:

- simulate the pieces separately (cavity with the qubit chip removed → bare cavity; isolated
  transmon → bare qubit), or
- de-hybridize: fit the two eigenmode frequencies to a 2×2 coupling model
  [[ω_q^bare, g], [g, ω_c^bare]] and back out ω^bare and g.

In the EPR philosophy ω_m^bare and g are not separately physical — only the hybridized eigenmodes
are — which is why they aren't produced directly.

## Summary

- ω_q, ω_c (paper) = **linear hybridized eigenmodes**: include g, exclude Lamb shift; intermediate.
- **Measured** qubit/resonator frequency = ω_m − Δ_m (add the Lamb shift). This ≠ ω_q, ω_c.
- **Bare uncoupled** (g = 0) = neither; requires a separate sim or de-hybridization.

**Read in the paper:** the "intermediate parameter" caveat after Eq. (16/17); the g "fully factored"
remark in Section II.A; the "dressed mode frequencies ω_m − Δ_m" in the Comparison section.
