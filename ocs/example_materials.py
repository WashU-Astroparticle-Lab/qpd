"""
Example demonstrating different superconducting materials in OCS

This script shows how to use the materials.yaml database to simulate
OCS transmons with different superconductors (Al, Hf, Nb, TiN)
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from ocs_transmon import OCS
from pathlib import Path

# Load custom style
_style_path = Path(__file__).parent / "ocs.mplstyle"


def example_compare_materials():
    """
    Compare OCS transmon behavior with different materials
    """
    print("\n" + "=" * 70)
    print("Material Comparison: OCS Transmon with Different Superconductors")
    print("=" * 70)
    
    # Same device parameters
    e_j_hz = 8.335e9  # ~8.3 GHz
    e_c_hz = 0.695e9  # ~695 MHz
    
    # Materials to compare
    materials = ['aluminum', 'hafnium', 'aluminum_manganese']
    
    # List all available materials
    print(f"\nAvailable materials: {OCS.list_materials()}")
    print(f"\nComparing: {materials}\n")
    
    # Create instances for each material
    ocs_instances = {}
    for mat in materials:
        ocs_instances[mat] = OCS(e_j_hz, e_c_hz, material=mat)
        props = ocs_instances[mat]
        print(f"{mat.upper()}:")
        print(f"  Tc = {props.tc:.3f} K")
        print(f"  Δ = {props.delta_material*1e6:.3f} μeV")
        print(f"  DOS = {props.dos:.2e} [1/(μm³·eV)]")
        print()
    
    # Compare energy level dispersion
    offset_charges = np.linspace(0, 1, 400)
    
    with plt.style.context(_style_path):
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7, 2.72))
        
        # Use GnBu colormap for odd parity comparisons
        cmap_odd = cm.get_cmap('GnBu')
        n_materials = len(materials)
        
        for idx, mat in enumerate(materials):
            ocs = ocs_instances[mat]
            energies_even, energies_odd, _ = ocs.solve_system(
                offset_charges, num_levels=2
            )
            
            # Get color from colormap
            color = cmap_odd(0.3 + 0.7 * idx / max(n_materials - 1, 1))
            
            # Plot f_01 for odd parity
            freq_01 = ((energies_odd[:, 1] - energies_odd[:, 0]) / 
                      ocs.PLANCK_EV_S / 1e9)
            
            # Calculate dispersion (max - min frequency)
            dispersion = np.max(freq_01) - np.min(freq_01)
            
            ax1.plot(offset_charges, freq_01, linewidth=2, color=color,
                    label=f'{mat} (Tc={ocs.tc:.2f}K)')
            
            # Plot ground state energy
            e_ground = (energies_odd[:, 0] - energies_odd[0, 0]) / (
                ocs.PLANCK_EV_S * 1e9
            )
            ax2.plot(offset_charges, e_ground * 1e3, linewidth=2,
                    color=color, label=f'{mat}')
        
        ax1.set_xlabel(r'Offset Charge [$C_g V_g / 2e$]')
        ax1.set_ylabel(r'$f_{01}$ [GHz] (odd parity)')
        ax1.set_title('Transition Frequency vs Offset Charge')
        ax1.legend()
        ax1.grid(alpha=0.3)
        
        ax2.set_xlabel(r'Offset Charge [$C_g V_g / 2e$]')
        ax2.set_ylabel(r'$E_0$ - $E_0$(u=0) [MHz]')
        ax2.set_title('Ground State Energy Dispersion')
        ax2.legend()
        ax2.grid(alpha=0.3)
    
    return fig


def example_custom_material():
    """
    Create OCS transmon with custom material properties
    """
    print("\n" + "=" * 70)
    print("Custom Material Properties Example")
    print("=" * 70)
    
    # Start with aluminum but override some properties
    e_j_hz = 8.335e9
    e_c_hz = 0.695e9
    
    print("\nDefault aluminum:")
    ocs_default = OCS(e_j_hz, e_c_hz, material='aluminum')
    print(f"  Tc = {ocs_default.tc} K")
    print(f"  DOS = {ocs_default.dos:.2e}")
    
    print("\nCustom aluminum (overridden Tc):")
    ocs_custom = OCS(e_j_hz, e_c_hz, material='aluminum', tc=0.8)
    print(f"  Tc = {ocs_custom.tc} K (overridden)")
    print(f"  DOS = {ocs_custom.dos:.2e} (from database)")
    
    return ocs_custom


def example_material_properties():
    """
    Display properties of all available materials
    """
    print("\n" + "=" * 70)
    print("Material Properties Database")
    print("=" * 70)
    
    materials = OCS.list_materials()
    
    print(f"\n{'Material':<20} {'Symbol':<8} {'Tc [K]':<10} "
          f"{'Δ [μeV]':<12} {'DOS [1/(μm³·eV)]'}")
    print("-" * 70)
    
    for mat in materials:
        props = OCS.get_material_properties(mat)
        symbol = mat[:2].capitalize() if mat != 'aluminum_manganese' else 'AlMn'
        
        # Convert to float in case YAML parsed as string
        tc = float(props['tc'])
        delta = float(props['delta'])
        dos = float(props['dos'])
        
        print(f"{mat:<20} {symbol:<8} {tc:<10.3f} "
              f"{delta*1e6:<12.3f} {dos:.2e}")
    
    print()


def example_hafnium_ocs():
    """
    Example using Hafnium for OCS transmon
    
    Hafnium has lower Tc than Al, which can be useful for studying
    quasiparticle effects at ultra-low temperatures.
    """
    print("\n" + "=" * 70)
    print("Hafnium OCS Transmon Example")
    print("=" * 70)
    
    # Hafnium device parameters
    e_j_hz = 6.0e9  # 6 GHz
    e_c_hz = 0.5e9  # 500 MHz
    
    ocs_hf = OCS(
        e_j_hz=e_j_hz,
        e_c_hz=e_c_hz,
        material='hafnium',
        temperature_k=0.015  # 15 mK
    )
    
    print(f"\nHafnium OCS transmon:")
    print(f"  Material: {ocs_hf.material_name}")
    print(f"  Tc = {ocs_hf.tc} K")
    print(f"  Δ = {ocs_hf.delta_material*1e6:.3f} μeV")
    print(f"  E_J/E_C = {ocs_hf.ej_ec_ratio:.1f}")
    print(f"  Temperature = {ocs_hf.temperature_k*1e3:.1f} mK")
    print(f"  T/Tc = {ocs_hf.temperature_k/ocs_hf.tc:.3f}")
    
    # Calculate some energy levels
    energies_even, energies_odd, energy_diff = ocs_hf.solve_system(
        [0, 0.25, 0.5], num_levels=3
    )
    
    print(f"\nEnergy levels (GHz):")
    for i, u in enumerate([0, 0.25, 0.5]):
        e01 = ((energies_odd[i, 1] - energies_odd[i, 0]) / 
              ocs_hf.PLANCK_EV_S / 1e9)
        print(f"  u={u:.2f}: E_01 = {e01:.3f} GHz")
    
    return ocs_hf


def main():
    """Run all material examples"""
    
    print("\n" + "=" * 70)
    print(" OCS Transmon - Material Properties Examples")
    print("=" * 70)
    
    # Display available materials
    example_material_properties()
    
    # Compare materials
    fig = example_compare_materials()
    
    # Custom properties
    ocs_custom = example_custom_material()
    
    # Hafnium example
    ocs_hf = example_hafnium_ocs()
    
    plt.show()
    
    print("\n" + "=" * 70)
    print("Material examples complete!")
    print("\nTo add new materials, edit materials.yaml")
    print("=" * 70 + "\n")


if __name__ == "__main__":
    main()

