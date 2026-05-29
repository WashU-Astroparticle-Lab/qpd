# Does the classical→quantum mapping (ω_m, φ_mj) assume tiny fields?

**Short answer:** Partly — but the assumption is *not* on the mapping itself.

- **Mapping is exact.** Extracting ω_m and φ_mj from the classical eigenmode solve is
  approximation-free. The only assumption is **linearizing each junction about its equilibrium**
  (junction → linear inductance L_j); ω_m, φ_mj are then the *exact* eigenmodes/ZPF of that
  linearized circuit. Eq. (21) is the exact bridge.
  - See: **Conclusion**, last paragraph ("derived within circuit theory without approximations…
    a change of basis adapted to nonlinear elements"); **Supplementary Section A**.

- **Where "tiny field" actually enters:** only when reducing Ĥ_nl to the analytic Kerr /
  cross-Kerr parameters — the **dispersive + weakly-anharmonic** regime,
  ω_k − ω_m ≫ E_J c_jp ⟨φ̂_j^p⟩ (p ≥ 3), absence of strong drives.
  - See: paragraph just before **Eq. (25)**.
  - The fully numerical diagonalization (pyEPR `cos_trunc`/`fock_trunc`) relaxes even this — it
    keeps the whole cosine.

**Takeaway:** the real assumption behind (ω_m, φ_mj) is *linearization about equilibrium*, not
small fields. Smallness is a condition for the *perturbative analytic formulas*, a separate step.
