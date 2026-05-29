# Why $|\psi_k\rangle \neq$ bare basis, is it a truncation artifact, and is $\chi$ ill-defined?

Follow-up to [09](09-chi-ND-labeling-bias.md). Three questions, taken in order.

## 1. Why is $|\psi_k\rangle$ different from the bare basis?

Bare basis = eigenstates of the *linear* part $\hat H_{\rm lin} = \sum_m \hbar\omega_m \hat a_m^\dagger\hat a_m$, i.e. product Fock states $|\vec n\rangle$. They diagonalize $\hat H_{\rm lin}$ but not $\hat H_{\rm ND} = \hat H_{\rm lin} + \hat H_{\rm nl}$. Split $\hat H_{\rm nl}$ by how it acts in the bare basis:

![diagonal vs off-diagonal split of H_nl](img/10_split.svg)

- **Diagonal (number-conserving) part** — all the Kerr / cross-Kerr terms — commutes with every $\hat n_m$. It leaves Fock states as **exact eigenstates** and only shifts their energies. *If $\hat H_{\rm nl}$ had only these, $|\psi_k\rangle = |\vec n\rangle$ exactly and $\chi_{\rm ND} = \chi_{\rm O1}$, with no labeling issue.*
- $|\psi_k\rangle \neq |\vec n\rangle$ comes **entirely from the off-diagonal, non-number-conserving terms** (counter-rotating pieces of the cosine). They connect different Fock states $\to$ eigenstates are superpositions $\to$ hybridization, avoided crossings, label swaps.

(Higher *number-conserving* Kerr like $\hat a^{\dagger 3}\hat a^3$ is still diagonal: it makes $\chi$ photon-number-dependent — doc 08's "which levels" issue — but does **not** hybridize. Only off-diagonal terms do.)

## 2. Is it different without truncation?

**Yes — not a truncation artifact.** The off-diagonal terms are genuinely in the exact $\cos\hat\varphi_j$ (all powers of $\hat a_m + \hat a_m^\dagger$). Truncation can add *spurious* extras if unconverged, but the core mixing is real. Cleanest proof: an **isolated** transmon's exact eigenstates are Mathieu functions, not harmonic Fock states — already dressed before any cavity or truncation. Couple a cavity $\to$ joint eigenstates are not product states. Exact, not numerical.

## 3. Is the dispersive shift ill-defined to begin with?

**No — but be precise about what is well-defined.**

- The **energy levels** $E(\psi_k)$ are exact and basis-independent. Always well-defined.
- The **dispersive shift** $\chi = E_{11} - E_{10} - E_{01} + E_{00}$ is well-defined whenever the four eigenstates map uniquely onto bare $|00\rangle$, $|10\rangle$, $|01\rangle$, $|11\rangle$ — i.e. the dressed$\leftrightarrow$bare labeling is an unambiguous bijection:

![chi well-defined iff unique labeling](img/10_welldef.svg)

That condition is *exactly the dispersive regime* — the regime the quantity is named for. There the eigenstates are dressed ($\neq$ bare) but still uniquely labelable, and $\chi$ is a robust spectral number.

**Mild dressing does NOT make $\chi$ ill-defined — $\chi$ is *supposed* to include the dressing.** The dispersive shift *is* a dressing effect; that the states aren't bare is the whole point.

It becomes **ambiguous** only when strong off-diagonal mixing (near a resonance) destroys the one-to-one labeling — when "the $|1,1\rangle$ state" no longer exists as a distinct dressed level. Then:

- the math isn't wrong — eigenstates and energies are still exact;
- but "the dispersive shift" is an **effective, reduced description that presupposes identifiable modes**, and that presupposition has failed. The shift goes photon-number- and mixing-dependent; no single $\chi$ captures it.

### Hierarchy

| Object | Status |
|---|---|
| eigenstates / spectrum | fundamental, always well-defined (exact) |
| dispersive shift $\chi$ | emergent parameter, well-defined **in its regime of validity** only |

Like "the frequency of mode m given n photons in mode n": meaningful and robust until the modes stop being separable, at which point you describe physics by the full spectrum, not by $\chi$.

This is also why `chi_O1` and `chi_ND` agree in the dispersive regime (off-diagonal terms perturbative, dressing weak, labeling clean) and diverge only where off-diagonal terms go resonant — i.e. where the concept itself frays.
