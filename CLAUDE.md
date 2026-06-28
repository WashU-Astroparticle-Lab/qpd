# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

QPD is a Python simulation library for **QPD Transmon** qubits, developed by the WashU Astroparticle Lab. QPD transmons operate in an intermediate regime (E_J/E_C ≈ 10–20) between Cooper-pair boxes and standard transmons, enabling direct measurement of charge parity through dispersive shifts. Key references: [Serniak et al., PRL 2019](https://arxiv.org/pdf/1903.00113), [arXiv:2405.17192](https://arxiv.org/pdf/2405.17192).

## Setup

```bash
pip install -e .
```

This installs the `qpd` package in editable mode. Core dependencies (numpy, scipy, matplotlib, pyyaml, iminuit, qutip) are declared in `pyproject.toml`. `qutip` is required by the driven-dissipative solver `qpd.theory.driven` (imported at module load), so it is a core dependency, not optional. The `resonator_tools` circle-fit (Probst notch) used by `notebooks/simulation.ipynb` §9 and the `daq` readout chain is not on PyPI; install it via `pip install ".[fit]"` (git source) or `pip install -e /path/to/resonator_tools`.

## Usage

```python
from qpd import QPD

qpd = QPD(e_j_hz=8.335e9, e_c_hz=0.695e9)
qpd.plot_all()
```

All three import forms work:
- `from qpd import QPD`
- `from qpd.theory import QPD`
- `from qpd.theory.transmon import QPD`

## Running Examples

```bash
python examples/example_usage.py
python examples/example_materials.py
```

There is no formal test suite or linter configured.

## Repository Structure

```
qpd/
├── pyproject.toml                     # Package metadata and dependencies
├── CLAUDE.md, README.md, LICENSE
├── src/qpd/                           # Installable package (src layout)
│   ├── __init__.py                    # Re-exports QPD
│   └── theory/                        # Theory/simulation subpackage
│       ├── __init__.py                # Re-exports QPD
│       ├── transmon.py                # Core module (~1177 lines, class QPD)
│       ├── materials.yaml             # Superconductor properties database
│       └── qpd.mplstyle              # PRL publication-quality plot style
├── examples/
│   ├── example_usage.py
│   └── example_materials.py
├── notebooks/
│   ├── qpd.ipynb
│   └── std_analysis.ipynb
└── docs/
    ├── theory.md                      # API, physics, and usage docs
    └── materials.md                   # Material system guide
```

## Architecture

All core code lives in `src/qpd/theory/transmon.py`, which defines the `QPD` class:

- **Construction**: `QPD(e_j_hz, e_c_hz, ...)` or factory `QPD.from_capacitance(...)`.
- **Core physics**: `build_hamiltonian()` → `solve_eigensystem()` (scipy.linalg.eigh) → `solve_system()` computes even/odd parity energy levels. `compute_dispersive_matrix()` calculates dispersive shifts χ and matrix elements.
- **Visualization**: `plot_energy_levels()`, `plot_matrix_elements()`, `plot_dispersive_shift()`, `plot_parity_shift_vs_frequency()`, `plot_parity_shift_vs_ng()`, `plot_all()`. All use `qpd.mplstyle` for PRL publication-quality figures. The canonical style path is `QPD._style_path`; all code (examples, notebooks) should reference this rather than constructing paths independently.
- **C_Q fitting**: `compute_quantum_capacitance()` is the forward primitive. `fit_quantum_capacitance()` runs an iminuit (MIGRAD + HESSE) least-squares fit with optional `fit_scale` / `fit_baseline` / `fit_offset` toggles; pass `profile_ratio=True` to additionally run `profile_ratio_likelihood()` and attach the asymmetric 1σ/2σ/3σ profile-likelihood interval on E_J/E_C under `fit['ratio_profile']`. `plot_capacitance_fit()` auto-picks the profile interval for the title when present (toggle via `use_profile`); `plot_likelihood_landscape()` plots the 2D joint (E_J, E_C) chi² contours.
- **Materials**: `materials.yaml` stores superconductor properties (Al, Hf, AlMn, Nb, TiN). Loaded on demand via `load_materials_database()`.

The `qpd.theory` subpackage is intended to house simulation/theory code; the repo may grow to include other subpackages (e.g., `qpd.fitting`) alongside it.

## Key Conventions

- Physical energies are in Hz (not eV or GHz) at the API boundary; internal computations convert as needed.
- The Cooper-pair-box Hamiltonian is periodic in `n_g` with period 1 (one Cooper pair = 2e). Any `offset_charges` array passed to `compute_quantum_capacitance` / `fit_quantum_capacitance` must already be normalized so that one charge period = Δn_g = 1. `cq_measured` units are arbitrary — the fitter's `scale` parameter absorbs the prefactor — but the x-axis scale is the one piece of calibration the fitter cannot infer for you.
- Vectorized NumPy operations with broadcasting are used throughout for array-based parameter scans.
- Jupyter notebooks in `notebooks/` are used for interactive analysis and are often large binary blobs in git.
- Documentation lives in `docs/`.
