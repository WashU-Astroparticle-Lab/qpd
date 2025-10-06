# Quick Start Guide

## Installation

```bash
cd OCS
pip install -r requirements.txt
```

**Note**: All plots use a custom style (`ocs.mplstyle`) for 
publication-quality figures with PRL-standard dimensions.

## Basic Usage (3 lines)

```python
from ocs_transmon import OCS
ocs = OCS(e_j_hz=8.335e9, e_c_hz=0.695e9)  # 8.3 GHz, 695 MHz
ocs.plot_all()
```

## Run Examples

```bash
python example_usage.py
```

Choose from menu:
1. WashU parameters
2. Serniak parameters  
3. From capacitance
4. Custom analysis
5. Parameter scan
6. Run all

## Verify Installation

```bash
python verify_translation.py
```

Should output transition frequencies, dispersive shifts, and verify 
symmetries.

## Common Recipes

### Recipe 1: Custom Parameters

```python
from ocs_transmon import OCS

# Your parameters (in Hz)
e_j_hz = 6.26e9  # ~6.3 GHz (equivalent to 0.3 K·kB)
e_c_hz = 0.417e9  # ~417 MHz (equivalent to 0.02 K·kB)

ocs = OCS(e_j_hz, e_c_hz, temperature_k=0.020, r_n_ohm=30e3)

# Single plot
fig, ax = ocs.plot_energy_levels()

# Or all plots
figs = ocs.plot_all(coupling_g_hz=200e6, readout_freq_hz=6.5e9)
```

### Recipe 2: From Capacitance

```python
# Get material properties
al_props = OCS.get_material_properties('aluminum')
delta_hz = al_props['delta'] / OCS.PLANCK_EV_S  # Convert eV to Hz

ocs = OCS.from_capacitance(
    c_total_f=70e-15,  # 70 fF
    delta_hz=delta_hz,  # Al gap in Hz
    r_n_ohm=25e3  # 25 kΩ
)
```

### Recipe 3: Compute Specific Values

```python
# Energy levels at offset charge u=0.25
energies_even, energies_odd, energy_diff = ocs.solve_system(
    offset_charges=[0.25], 
    num_levels=4
)

# Dispersive shift
matrix_elements, chi_ip = ocs.compute_dispersive_matrix(
    offset_charge=0.5,
    coupling_g_hz=150e6,
    readout_freq_hz=7e9,
    num_levels=6
)

print(f"χ₀ = {chi_ip[0] / 1e6:.3f} MHz")
```

### Recipe 4: Scan Parameter Space

```python
import numpy as np
import matplotlib.pyplot as plt

# Scan E_J/E_C ratios
ratios = [5, 10, 15, 20, 30, 50]
offset_charges = np.linspace(0, 1, 200)
e_j_hz = 8.335e9  # Fix E_J at ~8.3 GHz

for ratio in ratios:
    e_c_hz = e_j_hz / ratio
    ocs = OCS(e_j_hz, e_c_hz)
    
    _, energies_odd, _ = ocs.solve_system(offset_charges)
    freq_01 = (energies_odd[:, 1] - energies_odd[:, 0]) / (
        ocs.PLANCK_EV_S * 1e9
    )
    
    plt.plot(offset_charges, freq_01, label=f'{ratio}')

plt.xlabel('Offset Charge')
plt.ylabel('f₀₁ [GHz]')
plt.legend(title='$E_J/E_C$')
plt.show()
```

## Physical Constants

Access via class attributes:
```python
OCS.PLANCK_EV_S    # Planck constant [eV·s]
OCS.KB_EV_K        # Boltzmann constant [eV/K]
OCS.ELECTRON_CHARGE # Elementary charge [C]

# Material properties are loaded from materials.yaml
OCS.list_materials()  # List available materials
OCS.get_material_properties('aluminum')  # Get specific material
```

## Key Methods

| Method | Purpose |
|--------|---------|
| `solve_eigensystem()` | Diagonalize Hamiltonian |
| `solve_system()` | Get even/odd parity levels |
| `compute_dispersive_matrix()` | Calculate χ and matrix elements |
| `plot_energy_levels()` | Energy diagram |
| `plot_matrix_elements()` | Transition matrix elements |
| `plot_dispersive_shift()` | χ vs offset charge |
| `plot_parity_shift_vs_frequency()` | χ vs resonator freq |
| `plot_all()` | All four plots |

## Typical Parameter Ranges

| Parameter | Typical Value | OCS Regime |
|-----------|---------------|------------|
| E_J/E_C | 10-20 | ✓ |
| E_J | 5-15 GHz | |
| E_C | 0.3-1.5 GHz | |
| Temperature | 10-30 mK | |
| R_n | 10-50 kΩ | |
| Coupling g | 50-300 MHz | |
| Resonator | 5-8 GHz | |

## Output

Running `ocs.plot_all()` generates:

1. **Energy levels** vs offset charge (even/odd parity)
2. **Matrix elements** |⟨j|n̂|0⟩| (log scale)
3. **Dispersive shift** χᵢ vs offset charge
4. **Parity contrast** Δχ vs resonator frequency

Plus console output with key frequencies and parameters.

## Troubleshooting

**Import error?**
```bash
pip install numpy scipy matplotlib
```

**Values seem wrong?**
- Check units: energies in eV, frequencies in Hz
- Verify E_J/E_C in range 5-50 for OCS regime
- Try `verify_translation.py` for known good values

**Plots not showing?**
```python
import matplotlib.pyplot as plt
# ... your code ...
plt.show()  # Add this at the end
```

## Documentation

- Full API: See `README.md`
- Examples: See `example_usage.py`
- Translation: See `TRANSLATION_NOTES.md`
- Physics: See arxiv.org/pdf/1903.00113

## Support

For physics questions: Refer to Serniak et al. (2019)  
For code issues: Check `README.md` and `example_usage.py`

