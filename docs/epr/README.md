# EPR Study Notes

Notes on **Energy-Participation Quantization of Josephson Circuits**
([Minev et al., arXiv:2010.00620](https://arxiv.org/abs/2010.00620)), with cross-references
to the [`qiskit-metal`](https://github.com/qiskit-community/qiskit-metal) / `pyEPR` implementation.

Each file answers one question, focused on the physical picture. Equations are written in
KaTeX-compatible `$...$` / `$$...$$` so they render in Cursor's markdown preview.

## Index

_(entries added as questions are answered)_

| # | Topic | File |
|---|-------|------|
| 01 | Full Hamiltonian, term by term | [01-full-hamiltonian.md](01-full-hamiltonian.md) |
| 02 | Deriving Eq. (7) from Eq. (6): the ZPF identity φ²=p·ħω/2E_J | [02-eq6-to-eq7.md](02-eq6-to-eq7.md) |
| 03 | Does the ω_m, φ_mj mapping assume tiny fields? (linearization, not smallness) | [03-validity-tiny-field.md](03-validity-tiny-field.md) |
| 04 | Do Eq. (11) and Eq. (12) contradict? (no — single-junction Cauchy–Schwarz equality) | [04-eq11-eq12-consistency.md](04-eq11-eq12-consistency.md) |
| 05 | Seeing χ_qc as the dispersive shift from Eq. (8); Lamb vs dispersive; vacuum vs real | [05-chi-as-dispersive-shift.md](05-chi-as-dispersive-shift.md) |
| 06 | Does the AC Stark effect exist in Eq. (8)? (it's the cross-Kerr; no drive term) | [06-ac-stark-in-eq8.md](06-ac-stark-in-eq8.md) |
| 07 | Does χ_ND assume Eq. (8)? (no — ND is non-perturbative; its own limits) | [07-chi-ND-vs-perturbation.md](07-chi-ND-vs-perturbation.md) |
| 08 | Writing H_ND with α_m and χ_mn (Kerr form = Eq. 25; why it's only ≃) | [08-H_ND-in-kerr-form.md](08-H_ND-in-kerr-form.md) |
| 09 | χ_ND labeling: worked example + why it systematically biases dispersion | [09-chi-ND-labeling-bias.md](09-chi-ND-labeling-bias.md) |
| 10 | Why |ψ⟩≠bare basis; not a truncation artifact; is χ ill-defined? | [10-why-dressed-states-differ.md](10-why-dressed-states-differ.md) |
| 11 | Is the measured dispersive shift the same as χ_ND? (yes — same spectral observable) | [11-measured-vs-chiND.md](11-measured-vs-chiND.md) |
| 12 | Getting bare qubit/resonator frequencies; why they aren't ω_q, ω_c | [12-bare-vs-dressed-frequencies.md](12-bare-vs-dressed-frequencies.md) |
| 13 | **Summary:** extracting freq / anharmonicity / Lamb shift / dispersion from EPR (+pyEPR fields) | [13-extracting-parameters-summary.md](13-extracting-parameters-summary.md) |
| 14 | Intro to the Kerr Hamiltonian — demystifying self-Kerr/cross-Kerr/anharmonicity jargon | [14-kerr-hamiltonian-intro.md](14-kerr-hamiltonian-intro.md) |
| 15 | Physical origin of the Lamb shift (vacuum dressing via the nonlinearity) | [15-lamb-shift-physical-origin.md](15-lamb-shift-physical-origin.md) |

## Rendering note

Warp's markdown viewer has **no math engine**, so equations are pre-rendered to SVG images
(white card, legible on dark theme) under `img/`. Regenerate with `render_eqs.py`.
Doc 01 still uses raw `$...$` (works in Cursor/VS Code preview, not Warp); doc 02 onward uses images.
