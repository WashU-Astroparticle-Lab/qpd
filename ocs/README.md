# OCS Transmon Simulator

Python implementation for simulating **Offset-Charge-Sensitive (OCS) 
Transmons**, based on the physics described in:

- **Serniak et al., PRL (2019)**: [Direct Dispersive Monitoring of 
Charge Parity in Offset-Charge-Sensitive Transmons](
https://arxiv.org/pdf/1903.00113)
- Additional reference: [arXiv:2405.17192](
https://arxiv.org/pdf/2405.17192)

## Overview

OCS transmons operate in the intermediate regime between Cooper-pair 
box (E_J/E_C ≈ 1) and standard transmon (E_J/E_C ≳ 50), where charge 
dispersion is measurable. This makes them useful for:

- Probing charge parity dynamics
- Studying nonequilibrium quasiparticle (QP) effects
- Investigating superconducting gap variations
- Testing modified Josephson current-phase relations

## Installation

```bash
pip install -r requirements.txt
```

## Plotting Style

All plots use a custom matplotlib style (`ocs.mplstyle`) designed for 
publication-quality figures. The style features:
- Serif fonts (DejaVu Serif) with consistent sizing
- PRL-standard figure dimensions (3.375" × 2.72")
- Professional tick marks (inward-facing) and grid styling
- High-resolution output (600 DPI for saved figures)
- Constrained layout for optimal spacing

The style is automatically applied to all OCS plotting methods. To use 
it in custom scripts:

```python
from pathlib import Path
import matplotlib.pyplot as plt

style_path = Path("OCS") / "ocs.mplstyle"
with plt.style.context(style_path):
    # Your plotting code here
    fig, ax = plt.subplots()
    # ...
```

## Quick Start

```python
from ocs_transmon import OCS

# Initialize with Josephson and charging energies
kb = OCS.KB_EV_K  # Boltzmann constant in eV/K
e_j = 0.4 * kb  # Josephson energy
e_c = e_j / 12  # Charging energy (E_J/E_C = 12)

ocs = OCS(
    e_josephson=e_j,
    e_charging=e_c,
    temperature=0.02,  # 20 mK
    normal_resistance=27e3  # 27 kΩ
)

# Generate all standard plots
figs = ocs.plot_all(
    coupling_g=150e6,  # 150 MHz transmon-resonator coupling
    resonator_freq=7.0e9,  # 7 GHz resonator
    num_levels=6
)
```

## Features

### Core Functionality

The `OCS` class provides:

1. **Hamiltonian Construction**: Build Cooper-pair box Hamiltonian in 
charge basis
2. **Eigenspectrum Calculation**: Solve for energy levels with both 
even and odd charge parity
3. **Dispersive Shifts**: Compute charge-parity-dependent dispersive 
shifts χ_i
4. **Matrix Elements**: Calculate transition matrix elements for 
charge operator
5. **Visualization**: Generate publication-quality plots

### Available Methods

#### Initialization

```python
# Method 1: From E_J and E_C directly
ocs = OCS(e_josephson, e_charging, **kwargs)

# Method 2: From capacitance and resistance
ocs = OCS.from_capacitance(
    total_capacitance,  # C_Σ = C_J + C_g + C_shunt
    delta,  # Superconducting gap
    normal_resistance,
    **kwargs
)
```

#### Core Computation Methods

- `build_hamiltonian(offset_charge, charge_cutoff)`: Construct 
Hamiltonian matrix
- `solve_eigensystem(offset_charge, charge_cutoff)`: Diagonalize 
Hamiltonian
- `solve_system(offset_charges, num_levels)`: Compute even/odd parity 
levels
- `compute_dispersive_matrix(offset_charge, coupling_g, 
resonator_freq, num_levels)`: Calculate χ and matrix elements

#### Plotting Methods

- `plot_energy_levels()`: Energy diagram vs offset charge
- `plot_matrix_elements()`: Charge matrix elements |⟨j|n̂|0⟩|
- `plot_dispersive_shift()`: Dispersive shift χ vs offset charge
- `plot_parity_shift_vs_frequency()`: Parity-dependent χ vs frequency
- `plot_all()`: Generate all four standard plots

## Examples

See `example_usage.py` for comprehensive examples including:

1. **Basic usage** with WashU parameters
2. **Serniak parameters** from the reference paper
3. **Initialization from capacitance**
4. **Custom analysis workflows**
5. **Parameter scans** over E_J/E_C ratio

Run the examples:

```bash
python example_usage.py
```

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

## Key Parameters

### Physical Constants

Imported from `scipy.constants`:
- Planck's constant: `h`
- Boltzmann constant: `k`
- Elementary charge: `e`
- Electron-volt: `eV`

### Material Properties (Aluminum)

- Density of states: `1.72×10¹⁰ μm⁻³·eV⁻¹`
- Fermi level: `11.6 eV`
- Critical temperature: `1.2 K`
- Superconducting gap: `1.89×10⁻⁴ eV`

### Typical Experimental Values

- Temperature: `0.02 K` (20 mK)
- Normal resistance: `27 kΩ`
- E_J/E_C ratio: `10-20` (OCS regime)
- Transmon-resonator coupling: `100-200 MHz`
- Readout resonator: `6-8 GHz`

## Translation from MATLAB

This Python implementation is translated from the original MATLAB 
scripts with the following improvements:

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

### Code Structure Improvements

1. **Object-oriented design**: All functionality in `OCS` class
2. **No hardcoded parameters**: All values passed as arguments
3. **Physical constants imported**: Uses `scipy.constants`
4. **Vectorized operations**: NumPy broadcasting instead of loops where 
possible
5. **Modular plotting**: Separate methods for each plot type
6. **Type hints and documentation**: Comprehensive docstrings
7. **Consistent naming**: Snake_case with descriptive names

### File Mapping

| MATLAB File | Python Equivalent |
|-------------|-------------------|
| `OCS_simple.m` | `main()` in `ocs_transmon.py` |
| `eigensystem.m` | `OCS.solve_eigensystem()` |
| `solvesystem.m` | `OCS.solve_system()` |
| `computeEcEj.m` | `OCS.from_capacitance()` |
| `dispermatrix.m` | `OCS.compute_dispersive_matrix()` |

## Output

The code generates four plots:

1. **Energy Level Diagram**: Shows f_0j vs offset charge for even 
(solid) and odd (dashed) parity
2. **Matrix Elements**: |⟨j,o|n̂|0,o⟩| on log scale showing transition 
strengths
3. **Dispersive Shift**: χ_i vs offset charge for ground and first 
excited states
4. **Parity Contrast**: |Δχ_0| vs resonator frequency showing optimal 
readout frequency

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

## License

See the LICENSE file in the repository root.

## Contributing

Contributions are welcome! Please ensure:
- Code follows PEP 8 style guidelines
- Variable names use descriptive snake_case
- All methods have comprehensive docstrings
- Physical units are clearly documented
- Lines do not exceed 100 characters

## Contact

For questions about the physics, refer to the papers above. For code 
issues, please open an issue in the repository.

