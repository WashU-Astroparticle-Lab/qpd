# Material Properties System Guide

## Overview

The OCS transmon simulator now supports multiple superconducting materials through a flexible YAML-based configuration system. Material properties are no longer hardcoded, making it easy to:

- Switch between different superconductors (Al, Hf, Nb, TiN, etc.)
- Add custom materials by editing `materials.yaml`
- Override specific properties for custom experiments

## Available Materials

| Material | Symbol | Tc [K] | Δ [μeV] | DOS [1/(μm³·eV)] | Notes |
|----------|--------|--------|---------|------------------|-------|
| aluminum | Al | 1.200 | 189 | 1.72×10¹⁰ | Standard for transmons |
| hafnium | Hf | 0.128 | 22.5 | 2.50×10¹⁰ | Lower Tc, recent OCS studies |
| aluminum_manganese | AlMn | 0.100 | 15.7 | 1.72×10¹⁰ | Very low Tc |
| niobium | Nb | 9.200 | 1400 | 1.80×10¹⁰ | High Tc, resonators |
| titanium_nitride | TiN | 3.000 | 440 | 2.22×10¹⁰ | Tunable properties |

## Basic Usage

### Using Default Material (Aluminum)

```python
from ocs_transmon import OCS

# Aluminum is the default
ocs = OCS(e_j_hz=8.335e9, e_c_hz=0.695e9)
print(f"Material: {ocs.material_name}, Tc = {ocs.tc} K")
```

### Specifying a Different Material

```python
# Use hafnium instead
ocs_hf = OCS(
    e_j_hz=8.335e9,
    e_c_hz=0.695e9,
    material='hafnium'
)

print(f"Material: {ocs_hf.material_name}")
print(f"Tc = {ocs_hf.tc} K")
print(f"Δ = {ocs_hf.delta_material*1e6:.2f} μeV")
print(f"DOS = {ocs_hf.dos:.2e} [1/(μm³·eV)]")
```

### Listing Available Materials

```python
from ocs_transmon import OCS

# List all materials in database
materials = OCS.list_materials()
print(f"Available: {materials}")

# Get properties for a specific material
props = OCS.get_material_properties('hafnium')
print(f"Hafnium Tc: {props['tc']} K")
```

## Advanced Usage

### Overriding Material Properties

You can override specific properties while using a material from the database:

```python
# Start with aluminum but use a custom Tc
ocs = OCS(
    e_j_hz=8.335e9,
    e_c_hz=0.695e9,
    material='aluminum',
    tc=0.8  # Override Tc to 0.8 K
)

print(f"Custom Tc: {ocs.tc} K")
print(f"DOS from database: {ocs.dos:.2e}")
```

Override multiple properties:

```python
ocs = OCS(
    e_j_hz=8.335e9,
    e_c_hz=0.695e9,
    material='aluminum',
    tc=0.9,
    dos=2.0e10,  # Custom density of states
    delta=2.0e-4  # Custom gap [eV]
)
```

### Custom Superconducting Gaps

Set left and right gaps explicitly (for asymmetric junctions):

```python
# Convert Hz to eV for gaps if needed
delta_hz = 45.7e9  # 45.7 GHz

ocs = OCS(
    e_j_hz=8.335e9,
    e_c_hz=0.695e9,
    material='aluminum',
    delta_l_hz=delta_hz,  # Left gap
    delta_r_hz=delta_hz * 0.9  # Right gap 10% smaller
)
```

## Adding New Materials

Edit `materials.yaml` to add new materials:

```yaml
materials:
  your_material:
    name: "Your Material Name"
    symbol: "Ym"
    description: "Brief description"
    properties:
      dos: 2.0e10  # Density of states [1/(μm³·eV)]
      fermi_level: 10.0  # Fermi energy [eV]
      tc: 1.5  # Critical temperature [K]
      delta: 2.2e-4  # SC gap at T=0 [eV]
      notes: "Additional information"
```

**Important YAML formatting rules:**
- Use scientific notation without inline comments (e.g., `dos: 1.72e10`)
- The code automatically converts to float, so both `1.72e10` and `17200000000` work
- Comments go on separate lines or at the end after adequate spacing

## Material Properties Reference

### Property Definitions

- **dos**: Density of states at Fermi level [1/(μm³·eV)]
  - Affects quasiparticle density calculations
  - Higher DOS → more quasiparticles at given temperature

- **fermi_level**: Fermi energy [eV]
  - Sets the energy scale for electronic properties
  - Typically 5-12 eV for common superconductors

- **tc**: Critical temperature [K]
  - Temperature below which material is superconducting
  - Affects thermal quasiparticle population

- **delta**: Superconducting gap at T=0 [eV]
  - Energy gap in density of states
  - Related to Tc by BCS theory: Δ ≈ 1.764 × kB × Tc
  - Can deviate from BCS for strong-coupling superconductors

### BCS Relationship

For weak-coupling superconductors:

```python
Δ = 1.764 * kB * Tc

# In eV
kb_ev = 8.617e-5  # eV/K
delta_ev = 1.764 * kb_ev * tc

# Example for Al (Tc = 1.2 K)
delta_al = 1.764 * 8.617e-5 * 1.2 ≈ 1.82e-4 eV
```

## Material Selection Guide

### When to Use Each Material

**Aluminum (Al)**
- Default choice for transmon qubits
- Well-characterized properties
- Good for temperatures ≥ 20 mK
- Use when: Standard fabrication, well-established processes

**Hafnium (Hf)**
- Lower Tc than Al
- Higher density of states
- Use when: Ultra-low temperature experiments (< 20 mK), studying quasiparticle effects with lower thermal QP background

**Aluminum-Manganese (AlMn)**
- Very low Tc (100 mK)
- Use when: Experiments requiring minimal superconducting gap, special low-energy applications

**Niobium (Nb)**
- High Tc (9.2 K)
- Use when: Resonators, applications where high Tc reduces quasiparticle issues, robustness to temperature fluctuations

**Titanium Nitride (TiN)**
- Intermediate Tc, tunable
- High kinetic inductance
- Use when: High-Q resonators, kinetic inductance applications

## Example: Comparing Materials

```python
import numpy as np
import matplotlib.pyplot as plt
from ocs_transmon import OCS

# Same device parameters
e_j_hz = 8.335e9
e_c_hz = 0.695e9

# Compare materials
materials = ['aluminum', 'hafnium', 'niobium']
offset_charges = np.linspace(0, 1, 400)

for mat in materials:
    ocs = OCS(e_j_hz, e_c_hz, material=mat)
    _, energies_odd, _ = ocs.solve_system(offset_charges, num_levels=2)
    
    freq_01 = (energies_odd[:, 1] - energies_odd[:, 0]) / (
        ocs.PLANCK_EV_S * 1e9
    )
    
    plt.plot(offset_charges, freq_01, 
            label=f'{mat} (Tc={ocs.tc:.2f}K)')

plt.xlabel('Offset Charge')
plt.ylabel('f_01 [GHz]')
plt.legend()
plt.show()
```

## Practical Considerations

### Temperature Effects

The curly N parameter (quasiparticle density) depends on material:

```python
curly_n = DOS * sqrt(2π * Δ * kB * T)
```

Lower Tc materials (like Hf) have:
- Smaller gap Δ
- Potentially fewer thermal quasiparticles at same T/Tc
- Different coherence properties

### Fabrication

Different materials require different:
- Deposition techniques
- Junction formation methods
- Oxidation/barrier layers
- Annealing procedures

Consult fabrication literature for specific material requirements.

## References

- **Aluminum**: Serniak et al., PRA (2019) - https://arxiv.org/pdf/1903.00113
- **Hafnium**: Recent OCS studies - https://arxiv.org/pdf/2405.17192
- **General**: Material properties from standard superconductor databases

## Quick Reference Card

```python
# List materials
OCS.list_materials()

# Get properties
props = OCS.get_material_properties('hafnium')

# Create with material
ocs = OCS(e_j_hz, e_c_hz, material='hafnium')

# Override properties
ocs = OCS(e_j_hz, e_c_hz, material='aluminum', tc=0.8)

# Custom gaps
ocs = OCS(e_j_hz, e_c_hz, delta_l_hz=50e9, delta_r_hz=45e9)

# Access properties
print(ocs.material_name, ocs.tc, ocs.dos, ocs.delta_material)
```

## Troubleshooting

**Material not found error:**
```
ValueError: Material 'xxx' not found. Available: aluminum, hafnium, ...
```
→ Check spelling, use `OCS.list_materials()` to see available options

**YAML parsing error:**
→ Check `materials.yaml` syntax, ensure proper indentation and no inline comments with scientific notation

**Type errors:**
→ The code automatically converts to float, but verify YAML values are properly formatted numbers

## Contributing Materials

To contribute new material data:
1. Verify properties from peer-reviewed sources
2. Add to `materials.yaml` following the existing format
3. Include references in notes field
4. Test with `OCS.get_material_properties('your_material')`
5. Submit pull request with documentation

For questions about physical values, consult:
- BCS theory for gap calculations
- Density of states from band structure calculations
- Critical temperatures from experimental measurements

