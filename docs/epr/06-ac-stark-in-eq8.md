# Does the AC Stark effect exist in Eq. (8)?

**Reframe:** the AC Stark *mechanism* IS in Eq. (8) — it's literally the cross-Kerr term
−χ_qc n̂_q n̂_c. **AC Stark shift = dispersive shift = the same term**, named by which mode you call
the "field":

![one term, two names](img/06_same_term.svg)

- **Dispersive (readout):** cavity frequency depends on qubit state, ω_c − χ_qc n̂_q.
- **AC Stark:** qubit frequency depends on cavity photon number, ω_q − χ_qc n̂_c.

The qubit frequency with a populated cavity:

![AC Stark shift](img/06_acstark.svg)

## What you're right about

Eq. (8) has **no drive term** (no â + â† linear-in-field piece) — it is the bare, undriven
Hamiltonian. Nothing *imposes* a cavity population, so with the cavity in vacuum (n̄_c = 0) there is
no Stark shift, only the **Lamb shift** (the "1" in the 2n_c + 1 split — see doc 05).

## What would be wrong

To say the effect *can't occur* here. The susceptibility χ_qc — the shift per photon — is fully
present in Eq. (8). Add a drive, or prepare a coherent cavity state with n̄_c photons, and the qubit
shifts by χ_qc·n̄_c. **Eq. (8) supplies the susceptibility; the drive supplies the photons.**

## Summary

| | mechanism (χ_qc n̂_q n̂_c) | drive (â + â†) | net shift |
|---|---|---|---|
| In Eq. (8)? | **yes** | **no** | only Lamb shift until cavity is populated |
| AC Stark | this term | needed to set n̄_c | χ_qc n̄_c |
| Dispersive readout | this term | not needed (qubit state, not drive) | χ_qc per qubit excitation |

So AC Stark is "absent" only in the sense that no drive is written to activate it; the coupling
that produces it is the cross-Kerr term already in Eq. (8).

---

## Are the dispersive shift and the AC Stark shift the same effect?

**Yes** — same single term −χ_qc n̂_q n̂_c, two readings depending on which operator you hold fixed.
Because n̂_q n̂_c is symmetric, neither mode is privileged:

> ∂ω_c/∂n_q = ∂ω_q/∂n_c = −χ_qc.

- **Dispersive shift** (resonator view): fix qubit state → cavity at ω_c − χ_qc n̂_q. Driven by the
  qubit's *own state*; no extra drive needed. Used for **readout**.
- **AC Stark shift** (qubit view): fix cavity photons → qubit at ω_q − χ_qc n̂_c. Needs *real
  cavity photons* (drive or thermal). Used for **spectroscopy / photon-number calibration**.

Two caveats so "same effect" doesn't over-merge them:

1. **What acts as "the field" differs.** Dispersive readout uses the qubit state; AC Stark needs
   real photons in the cavity. (This is why Eq. (8) has the coupling but no drive — see above.)
2. **Photon-number splitting** is the same coupling in the resolved limit χ_qc > κ (cavity
   linewidth): the qubit line splits into discrete peaks at ω_q − χ_qc·n_c, one per Fock state —
   the AC Stark shift "quantized." Still the same χ_qc.

**One coupling, one number χ_qc — dispersive shift, AC Stark shift, and photon-number splitting are
its different experimental faces.**
