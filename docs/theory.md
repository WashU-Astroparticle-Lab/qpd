# QPD Transmon Simulator

Python implementation for simulating **QPD Transmons**, based on the
physics described in:

- **Serniak et al., PRL (2019)**: [Direct Dispersive Monitoring of
Charge Parity in Offset-Charge-Sensitive Transmons](
https://arxiv.org/pdf/1903.00113)
- Additional reference: [arXiv:2405.17192](
https://arxiv.org/pdf/2405.17192)

## Overview

QPD transmons operate in the intermediate regime between Cooper-pair
box (E_J/E_C ≈ 1) and standard transmon (E_J/E_C ≳ 50), where charge
dispersion is measurable. This makes them useful for:

- Probing charge parity dynamics
- Studying nonequilibrium quasiparticle (QP) effects
- Investigating superconducting gap variations
- Testing modified Josephson current-phase relations

## Installation

```bash
pip install -r requirements.txt   # from repo root
```

## Quick Start

```python
from qpd import QPD

qpd = QPD(e_j_hz=8.335e9, e_c_hz=0.695e9)  # 8.3 GHz, 695 MHz
qpd.plot_all()
```

Run the interactive examples:

```bash
python examples/example_usage.py
```

## Plotting Style

All plots use a custom matplotlib style (`qpd.mplstyle`) designed for
publication-quality figures:
- Serif fonts (DejaVu Serif) with consistent sizing
- PRL-standard figure dimensions (3.375" × 2.72")
- Professional tick marks (inward-facing) and grid styling
- High-resolution output (600 DPI for saved figures)
- Constrained layout for optimal spacing

The style is automatically applied by all `QPD` plotting methods. To
use it in custom scripts, reference the canonical path via the class:

```python
import matplotlib.pyplot as plt
from qpd import QPD

with plt.style.context(QPD._style_path):
    fig, ax = plt.subplots()
    # Your plotting code here
```

## API Reference

### Initialization

```python
# Method 1: From E_J and E_C directly (in Hz)
qpd = QPD(e_j_hz, e_c_hz, temperature_k=0.02, r_n_ohm=27e3)

# Method 2: From capacitance and resistance
qpd = QPD.from_capacitance(
    c_total_f=70e-15,   # C_Σ = C_J + C_g + C_shunt [F]
    delta_hz=delta_hz,   # Superconducting gap [Hz]
    r_n_ohm=25e3
)
```

### Core Computation Methods

| Method | Purpose |
|--------|---------|
| `build_hamiltonian(offset_charge, charge_cutoff)` | Construct CPB Hamiltonian matrix |
| `solve_eigensystem(offset_charge, charge_cutoff)` | Diagonalize Hamiltonian |
| `solve_system(offset_charges, num_levels)` | Compute even/odd parity energy levels |
| `compute_dispersive_matrix(offset_charge, coupling_g_hz, readout_freq_hz, num_levels)` | Calculate χ and matrix elements |

### Plotting Methods

| Method | Purpose |
|--------|---------|
| `plot_energy_levels()` | Energy diagram vs offset charge |
| `plot_matrix_elements()` | Charge matrix elements \|⟨j\|n̂\|0⟩\| |
| `plot_dispersive_shift()` | Dispersive shift χ vs offset charge |
| `plot_parity_shift_vs_frequency()` | Parity-dependent χ vs frequency |
| `plot_parity_shift_vs_ng()` | Parity shift vs gate charge |
| `plot_all()` | Generate all standard plots |

### Physical Constants & Materials

```python
QPD.PLANCK_EV_S       # Planck constant [eV·s]
QPD.KB_EV_K           # Boltzmann constant [eV/K]
QPD.ELECTRON_CHARGE   # Elementary charge [C]

QPD.list_materials()                        # Available materials
QPD.get_material_properties('aluminum')     # Get specific material
```

## Common Recipes

### From Capacitance

```python
al_props = QPD.get_material_properties('aluminum')
delta_hz = al_props['delta'] / QPD.PLANCK_EV_S

qpd = QPD.from_capacitance(
    c_total_f=70e-15,
    delta_hz=delta_hz,
    r_n_ohm=25e3
)
```

### Compute Specific Values

```python
# Energy levels at offset charge u=0.25
energies_even, energies_odd, energy_diff = qpd.solve_system(
    offset_charges=[0.25],
    num_levels=4
)

# Dispersive shift
matrix_elements, chi_ip = qpd.compute_dispersive_matrix(
    offset_charge=0.5,
    coupling_g_hz=150e6,
    readout_freq_hz=7e9,
    num_levels=6
)

print(f"χ₀ = {chi_ip[0] / 1e6:.3f} MHz")
```

### Scan Parameter Space

```python
import numpy as np

ratios = [5, 10, 15, 20, 30, 50]
offset_charges = np.linspace(0, 1, 200)
e_j_hz = 8.335e9

for ratio in ratios:
    e_c_hz = e_j_hz / ratio
    qpd = QPD(e_j_hz, e_c_hz)
    _, energies_odd, _ = qpd.solve_system(offset_charges)
    freq_01 = (energies_odd[:, 1] - energies_odd[:, 0]) / (
        qpd.PLANCK_EV_S * 1e9
    )
    plt.plot(offset_charges, freq_01, label=f'{ratio}')
```

## Typical Parameter Ranges

| Parameter | Typical Value | QPD Regime |
|-----------|---------------|------------|
| E_J/E_C | 10-20 | ✓ |
| E_J | 5-15 GHz | |
| E_C | 0.3-1.5 GHz | |
| Temperature | 10-30 mK | |
| R_n | 10-50 kΩ | |
| Coupling g | 50-300 MHz | |
| Resonator | 5-8 GHz | |

## Physics Background

### Cooper-Pair Box Hamiltonian

The system is described by:

```
H = 4E_C(n̂ - n_g)² - (E_J/2)(|n⟩⟨n+1| + h.c.)
```

where:
- `E_C = e²/(2C_Σ)`: Charging energy
- `E_J`: Josephson energy
- `n̂`: Cooper pair number operator
- `n_g = C_g V_g / (2e)`: Dimensionless offset charge

### Charge Parity

- **Even parity**: Integer number of Cooper pairs (ground state energy
at n_g)
- **Odd parity**: Unpaired electron present (ground state energy at
n_g + 0.5)

The energy splitting between parities provides sensitivity to
quasiparticle tunneling.

### Dispersive Shift

The dispersive shift of a coupled resonator is:

```
χ_i = g² ∑_{j≠i} [2ω_{ij} |⟨j,p|n̂|i,p⟩|² / (ω_{ij}² - ω_r²)]
```

where:
- `g`: Transmon-resonator coupling
- `ω_{ij}`: Transition frequency between levels i and j
- `ω_r`: Resonator frequency
- `p`: Charge parity (even or odd)

The parity-dependent shift enables direct dispersive readout of charge
parity.

## Translation from MATLAB

### Variable Naming Conventions

| MATLAB | Python | Description |
|--------|--------|-------------|
| `Ec`, `Ej` | `e_c`, `e_j` | Charging/Josephson energy |
| `u` | `offset_charge(s)` | Offset charge |
| `EE`, `EO` | `energies_even`, `energies_odd` | Parity-resolved energies |
| `DE` | `energy_diff` | Parity splitting |
| `nlevels` | `num_levels` | Number of energy levels |
| `chi_ip` | `chi_ip` | Dispersive shift |
| `matrixelems` | `matrix_elements` | Transition matrix elements |
| `g` | `coupling_g` | Coupling strength |
| `fr`, `wr` | `resonator_freq` | Resonator frequency |
| `Rn` | `r_n` or `normal_resistance` | Normal resistance |
| `planck` | `PLANCK_EV_S` | Planck constant |
| `kb` | `KB_EV_K` | Boltzmann constant |

### File Mapping

| MATLAB File | Python Equivalent |
|-------------|-------------------|
| `OCS_simple.m` | `main()` in `qpd.theory.transmon.py` |
| `eigensystem.m` | `QPD.solve_eigensystem()` |
| `solvesystem.m` | `QPD.solve_system()` |
| `computeEcEj.m` | `QPD.from_capacitance()` |
| `dispermatrix.m` | `QPD.compute_dispersive_matrix()` |

## Quantum Capacitance and Readout-Basis Mapping

The package provides primitives for working with the quantum
capacitance `C_Q(n_g)` and for inverting resonator-basis readout
measurements.

### Quantum capacitance

For each parity branch, the quantum capacitance is the curvature of
the level energy with respect to gate charge:

```
C_Q[F] = -(C_g/2e)^2 · h · ∂²E_level[Hz] / ∂n_g²
```

where `n_g` is the dimensionless offset charge in units of 2e and
`E_level[Hz]` is the energy of the chosen level. The package computes
this numerically (no asymptotic approximation) from the same
Hamiltonian used elsewhere:

```python
cq_even, cq_odd = qpd.compute_quantum_capacitance(
    offset_charges, c_g_f=None, level=0,
)
```

If `c_g_f` is supplied, results are returned in Farads. If left as
`None`, the function returns the *intrinsic* curvature
`-∂²E[Hz]/∂n_g²` in Hz; this is the form the fitter consumes
internally and lets callers apply any preferred prefactor.

### Fitting a measured C_Q trace

`QPD.fit_quantum_capacitance(offset_charges, cq_measured, ...)` fits
the numerical model to a measured trace and extracts `(E_J, E_C, n_g0,
scale)`. The optional `n_g0` accounts for an imperfectly-centred
charge-bias axis; the optional `scale` absorbs any unknown
amplitude prefactor (so a known `C_g` is *not* required). When the
data is in a known absolute convention, pass `fit_scale=False,
fixed_scale=...` to break the (E_C, scale) shape/amplitude
degeneracy and get tight bounds on `E_J` and `E_C` individually.

### Phase shift → quantum capacitance (resonator basis)

In the resonator basis the qubit acts as a small capacitive load on
the cavity, so a measured frequency shift maps linearly to a change in
quantum capacitance:

```
Δf_r / f_r ≈ - ΔC_Q / (2 C_r)        (low-loading limit)
```

A driven Lorentzian cavity has phase response
`φ(Δ_d) = arctan(2 Δ_d / κ)` where `Δ_d = f_drive − f_r`. A *change*
in measured phase Δφ then corresponds to a frequency shift

```
Δf_r = Δ_d − (κ/2) · tan(arctan(2 Δ_d / κ) + Δφ)
```

The two pure helpers in `qpd.theory.readout` implement both steps and
work for arbitrary drive placement (the readout drive does not have to
sit on either parity-resolved resonance):

```python
from qpd import (
    phase_shift_to_frequency_shift,
    frequency_shift_to_quantum_capacitance,
    kappa_from_quality_factors,
)

kappa   = kappa_from_quality_factors(f_r, q_i, q_c)   # Q_i, Q_c -> κ
df_r    = phase_shift_to_frequency_shift(delta_phi, kappa, drive_detuning)
delta_q = frequency_shift_to_quantum_capacitance(df_r, f_r, c_r)
```

The κ helper uses the standard parallel-Q rule
`1/Q_loaded = 1/Q_i + 1/Q_c`, `κ = f_r / Q_loaded`. Convention:
real `Q_c` (symmetric hanger / notch). For asymmetric hangers
where Q_c is complex, compute the loaded Q yourself and pass
`κ = f_r / Q_loaded` directly.

See `examples/example_capacitance.py` for an end-to-end demonstration
covering forward computation, the round-trip fit, and the readout
chain.

## References

1. K. Serniak et al., "Direct Dispersive Monitoring of Charge Parity
in Offset-Charge-Sensitive Transmons," Phys. Rev. Lett. (2019).
[arXiv:1903.00113](https://arxiv.org/pdf/1903.00113)

2. Additional reference: [arXiv:2405.17192](
https://arxiv.org/pdf/2405.17192)

3. J. Koch et al., "Charge-insensitive qubit design derived from the
Cooper pair box," Phys. Rev. A 76, 042319 (2007).

4. G. Catelani, "Parity switching and decoherence by quasiparticles in
single-junction transmons," Phys. Rev. B 89, 094522 (2014).
