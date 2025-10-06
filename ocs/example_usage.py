"""
Example usage of the OCS transmon simulator

This script demonstrates various ways to use the OCS class to analyze
offset-charge-sensitive transmons.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from ocs_transmon import OCS
from pathlib import Path

# Load custom style
_style_path = Path(__file__).parent / "ocs.mplstyle"


def example_1_basic_usage():
    """
    Example 1: Basic usage with specified Eⱼ and Eᴄ
    
    This replicates the main MATLAB script behavior.
    """
    print("\n" + "=" * 70)
    print("Example 1: Basic OCS Transmon Simulation (WashU parameters)")
    print("=" * 70)
    
    # WashU parameters from MATLAB script (E_J/E_C = 15)
    ej_ec_ratio = 15
    e_j_hz = 8.335e9  # ~8.3 GHz (equivalent to 0.4 K·kB)
    e_c_hz = e_j_hz / ej_ec_ratio  # Hz
    
    ocs = OCS(
        e_j_hz=e_j_hz,
        e_c_hz=e_c_hz,
        temperature_k=0.02,  # 20 mK
        r_n_ohm=16.5e3  # 16.5 kΩ
    )
    
    # Generate all plots
    figs = ocs.plot_all(
        coupling_g_hz=150e6,  # 150 MHz
        readout_freq_hz=7.0e9,  # 7 GHz
        num_levels=5
    )
    
    return ocs, figs


def example_2_serniak_parameters():
    """
    Example 2: Serniak et al. parameters from arXiv:1903.00113
    """
    print("\n" + "=" * 70)
    print("Example 2: Serniak Parameters (arXiv:1903.00113)")
    print("=" * 70)
    
    # Serniak parameters (commented in MATLAB script)
    # E_J = 0.295 K·kB ≈ 6.15 GHz, E_C = 0.017 K·kB ≈ 0.35 GHz
    e_j_hz = 6.149e9  # ~6.15 GHz (equivalent to 0.295 K·kB)
    e_c_hz = 0.354e9  # ~354 MHz (equivalent to 0.017 K·kB)
    
    ocs = OCS(
        e_j_hz=e_j_hz,
        e_c_hz=e_c_hz,
        temperature_k=0.02,
        r_n_ohm=27e3
    )
    
    # Generate all plots
    figs = ocs.plot_all(
        coupling_g_hz=150e6,
        readout_freq_hz=7.0e9,
        num_levels=5
    )
    
    return ocs, figs


def example_3_from_capacitance():
    """
    Example 3: Initialize from capacitance and resistance
    """
    print("\n" + "=" * 70)
    print("Example 3: Initialize from Capacitance")
    print("=" * 70)
    
    # Capacitances from MATLAB script
    c_j = 4.43e-16  # F, junction capacitance
    c_g = 1e-15  # F, gate capacitance
    c_shunt = 6.5e-14  # F, shunt capacitance
    c_total = c_j + c_g + c_shunt
    
    # Superconducting gap (convert from eV to Hz)
    delta_hz = OCS.DELTA_AL / OCS.PLANCK_EV_S  # Hz
    r_n = 27e3  # Normal resistance
    
    ocs = OCS.from_capacitance(
        c_total_f=c_total,
        delta_hz=delta_hz,
        r_n_ohm=r_n,
        temperature_k=0.02
    )
    
    print(f"\nComputed from capacitances:")
    print(f"Cᵨ = {c_total * 1e15:.3f} fF")
    print(f"Eⱼ = {ocs.e_j_hz / 1e9:.3f} GHz")
    print(f"Eᴄ = {ocs.e_c_hz / 1e9:.4f} GHz")
    print(f"Eⱼ/Eᴄ = {ocs.ej_ec_ratio:.2f}")
    
    return ocs


def example_4_custom_analysis():
    """
    Example 4: Custom analysis without using plot_all
    
    Shows how to access individual methods for custom workflows.
    """
    print("\n" + "=" * 70)
    print("Example 4: Custom Analysis Workflow")
    print("=" * 70)
    
    e_j_hz = 8.335e9  # ~8.3 GHz
    e_c_hz = e_j_hz / 12
    
    ocs = OCS(e_j_hz, e_c_hz, temperature_k=0.02, r_n_ohm=27e3)
    
    # Compute energy levels at specific offset charges
    offset_charges = np.array([0, 0.25, 0.5])
    energies_even, energies_odd, energy_diff = ocs.solve_system(
        offset_charges, num_levels=4
    )
    
    print("\nEnergy levels (GHz) at different offset charges:")
    print(f"{'Offset':<10} {'E₀₀(even)':<15} {'E₀₀(odd)':<15} "
          f"{'ΔE':<15}")
    print("-" * 60)
    for i, u in enumerate(offset_charges):
        e_even = energies_even[i, 0] / ocs.PLANCK_EV_S / 1e9
        e_odd = energies_odd[i, 0] / ocs.PLANCK_EV_S / 1e9
        delta_e = energy_diff[i] / ocs.PLANCK_EV_S / 1e9
        print(f"{u:<10.2f} {e_even:<15.6f} {e_odd:<15.6f} "
              f"{delta_e:<15.6f}")
    
    # Compute dispersive shift at a specific point
    coupling_g_hz = 150e6
    readout_freq_hz = 7.0e9
    matrix_elements, chi_ip = ocs.compute_dispersive_matrix(
        offset_charge=0.5,
        coupling_g_hz=coupling_g_hz,
        readout_freq_hz=readout_freq_hz,
        num_levels=4
    )
    
    print(f"\nDispersive shifts at offset charge 0.5:")
    print(f"χ₀ = {chi_ip[0] / 1e6:.3f} MHz")
    print(f"χ₁ = {chi_ip[1] / 1e6:.3f} MHz")
    
    print(f"\nMatrix elements |⟨j|n̂|0⟩|:")
    for j in range(1, 4):
        print(f"|⟨{j}|n̂|0⟩| = {np.sqrt(matrix_elements[0, j]):.6f}")
    
    return ocs


def example_5_parameter_scan():
    """
    Example 5: Scan over different Eⱼ/Eᴄ ratios
    
    Shows how the charge dispersion changes with Eⱼ/Eᴄ ratio.
    """
    print("\n" + "=" * 70)
    print("Example 5: Parameter Scan over Eⱼ/Eᴄ")
    print("=" * 70)
    
    e_j_hz = 8.335e9  # ~8.3 GHz
    
    # Scan different ratios
    ej_ec_ratios = [5, 10, 15, 20, 30, 50]
    offset_charges = np.linspace(0, 1, 200)
    
    with plt.style.context(_style_path):
        fig, ax = plt.subplots(figsize=(3.375, 2.72))
        
        # Use GnBu colormap for odd parity
        cmap_odd = cm.get_cmap('GnBu')
        n_ratios = len(ej_ec_ratios)
        
        for idx, ratio in enumerate(ej_ec_ratios):
            e_c_hz = e_j_hz / ratio
            ocs = OCS(e_j_hz, e_c_hz, temperature_k=0.02, r_n_ohm=27e3)
            
            energies_even, energies_odd, _ = ocs.solve_system(
                offset_charges, num_levels=2
            )
            
            # Get color from colormap
            color = cmap_odd(0.3 + 0.7 * idx / max(n_ratios - 1, 1))
            
            # Plot f_01 for odd parity
            freq_01 = ((energies_odd[:, 1] - energies_odd[:, 0]) / 
                      ocs.PLANCK_EV_S / 1e9)
            ax.plot(offset_charges, freq_01, linewidth=2, color=color,
                   label=f'$E_J/E_C = {ratio}$')
        
        ax.set_xlabel(r'Offset Charge [$C_g V_g / 2e$]')
        ax.set_ylabel(r'$f_{01}$ [GHz] (odd parity)')
        ax.set_title('Charge Dispersion vs $E_J/E_C$ Ratio')
        ax.legend()
        ax.grid(alpha=0.3)
        ax.minorticks_on()
    
    print(f"\nScanned Eⱼ/Eᴄ ratios: {ej_ec_ratios}")
    print("Higher ratios → smaller charge dispersion (more transmon-like)")
    print("Lower ratios → larger charge dispersion (more CPB-like)")
    
    return fig


def main():
    """Run all examples"""
    
    print("\n" + "=" * 70)
    print(" OCS Transmon Simulator - Example Usage")
    print("=" * 70)
    print("\nThis script demonstrates various uses of the OCS class.")
    print("Choose an example to run, or run all:\n")
    print("  1. Basic usage (WashU parameters)")
    print("  2. Serniak parameters")
    print("  3. Initialize from capacitance")
    print("  4. Custom analysis workflow")
    print("  5. Parameter scan over Eⱼ/Eᴄ")
    print("  6. Run all examples")
    
    choice = input("\nEnter choice (1-6) or press Enter for all: ")
    
    if choice == "1":
        ocs, figs = example_1_basic_usage()
    elif choice == "2":
        ocs, figs = example_2_serniak_parameters()
    elif choice == "3":
        ocs = example_3_from_capacitance()
    elif choice == "4":
        ocs = example_4_custom_analysis()
    elif choice == "5":
        fig = example_5_parameter_scan()
    else:
        # Run all examples
        print("\nRunning all examples...\n")
        ocs1, figs1 = example_1_basic_usage()
        ocs2, figs2 = example_2_serniak_parameters()
        ocs3 = example_3_from_capacitance()
        ocs4 = example_4_custom_analysis()
        fig5 = example_5_parameter_scan()
    
    plt.show()
    
    print("\n" + "=" * 70)
    print("Examples complete!")
    print("=" * 70 + "\n")


if __name__ == "__main__":
    main()

