# Does the classical$\to$quantum mapping ($\omega_m$, $\varphi_{mj}$) assume tiny fields?

**Short answer:** Partly — but the assumption is *not* on the mapping itself.

- **Mapping is exact.** Extracting $\omega_m$ and $\varphi_{mj}$ from the classical eigenmode solve is approximation-free. The only assumption is **linearizing each junction about its equilibrium** (junction $\to$ linear inductance $L_j$); $\omega_m$, $\varphi_{mj}$ are then the *exact* eigenmodes/ZPF of that linearized circuit. Eq. (21) is the exact bridge.
  - See: **Conclusion**, last paragraph ("derived within circuit theory without approximations… a change of basis adapted to nonlinear elements"); **Supplementary Section A**.

- **Where "tiny field" actually enters:** only when reducing $\hat H_{\rm nl}$ to the analytic Kerr / cross-Kerr parameters — the **dispersive + weakly-anharmonic** regime, $\omega_k - \omega_m \gg E_J c_{jp} \langle\hat\varphi_j^p\rangle$ ($p \ge 3$), absence of strong drives.
  - See: paragraph just before **Eq. (25)**.
  - The fully numerical diagonalization (pyEPR `cos_trunc`/`fock_trunc`) relaxes even this — it keeps the whole cosine.

**Takeaway:** the real assumption behind ($\omega_m$, $\varphi_{mj}$) is *linearization about equilibrium*, not small fields. Smallness is a condition for the *perturbative analytic formulas*, a separate step.
