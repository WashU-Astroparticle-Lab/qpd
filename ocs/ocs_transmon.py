"""
Offset-Charge-Sensitive (OCS) Transmon Analysis

This module implements simulation and analysis tools for OCS transmons,
based on the physics described in:
- Serniak et al., PRA (2019): https://arxiv.org/pdf/1903.00113
- Additional context from https://arxiv.org/pdf/2405.17192

The OCS transmon is in the intermediate regime between Cooper-pair box 
(E_J/E_C ≈ 1) and standard transmon (E_J/E_C ≳ 50), where charge 
dispersion is measurable.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.constants import h, k, e, eV
from scipy.linalg import eigh
import yaml
from pathlib import Path


class OCS:
    """
    Offset-Charge-Sensitive Transmon Simulator
    
    This class simulates OCS transmons by diagonalizing the Cooper-pair
    box Hamiltonian and computing dispersive shifts for charge-parity
    readout.
    
    Material properties are loaded from materials.yaml, supporting
    multiple superconductors (Al, Hf, Nb, TiN, etc.)
    """
    
    # Physical constants (in SI unless specified)
    PLANCK_EV_S = h / eV  # Planck constant in eV·s
    KB_EV_K = k / eV  # Boltzmann constant in eV/K
    ELECTRON_CHARGE = e  # Elementary charge in Coulombs
    
    # Materials database cache
    _materials_db = None
    _materials_path = Path(__file__).parent / "materials.yaml"
    _style_path = Path(__file__).parent / "ocs.mplstyle"
    
    @classmethod
    def load_materials_database(cls):
        """
        Load materials database from YAML file
        
        Returns
        -------
        dict
            Dictionary containing material properties
        """
        if cls._materials_db is None:
            with open(cls._materials_path, 'r') as f:
                cls._materials_db = yaml.safe_load(f)
        return cls._materials_db
    
    @classmethod
    def list_materials(cls):
        """
        List available materials in the database
        
        Returns
        -------
        list
            List of available material names
        """
        db = cls.load_materials_database()
        return list(db['materials'].keys())
    
    @classmethod
    def get_material_properties(cls, material_name):
        """
        Get properties for a specific material
        
        Parameters
        ----------
        material_name : str
            Name of the material (e.g., 'aluminum', 'hafnium')
            
        Returns
        -------
        dict
            Material properties dictionary
        """
        db = cls.load_materials_database()
        if material_name not in db['materials']:
            available = ', '.join(cls.list_materials())
            raise ValueError(
                f"Material '{material_name}' not found. "
                f"Available: {available}"
            )
        return db['materials'][material_name]['properties']
    
    def __init__(self, 
                 e_j_hz, e_c_hz, temperature_k=0.02, r_n_ohm=27e3, 
                 delta_l_hz=None, delta_r_hz=None,
                 material='aluminum', **material_overrides):
        """
        Initialize OCS transmon simulator
        
        Parameters
        ----------
        e_j_hz : float
            Josephson energy [Hz]
        e_c_hz : float
            Charging energy [Hz]
        temperature_k : float, optional
            Operating temperature [K], default 0.02 K
        r_n_ohm : float, optional
            Normal state resistance [Ω], default 27 kΩ
        delta_l_hz : float, optional
            Left superconducting gap [Hz], defaults to material value
        delta_r_hz : float, optional
            Right superconducting gap [Hz], defaults to material value
        material : str, optional
            Material name from database, default 'aluminum'
            Options: 'aluminum', 'hafnium', 'niobium', etc.
        **material_overrides : dict
            Override material properties (dos, fermi_level, tc, delta)
            Values should be in SI/eV units as specified in YAML
        """
        # Store input parameters in their native units
        self.e_j_hz = e_j_hz
        self.e_c_hz = e_c_hz
        self.temperature_k = temperature_k
        self.r_n_ohm = r_n_ohm
        
        # Load material properties
        self.material_name = material
        mat_props = self.get_material_properties(material).copy()
        
        # Apply any manual overrides
        mat_props.update(material_overrides)
        
        # Store material properties as instance variables
        # Convert to float in case YAML parsed as string
        self.dos = float(mat_props['dos'])  # [1/(μm³·eV)]
        self.fermi_level = float(mat_props['fermi_level'])  # [eV]
        self.tc = float(mat_props['tc'])  # [K]
        self.delta_material = float(mat_props['delta'])  # [eV]
        
        # Convert to eV for internal calculations
        self.e_j_ev = e_j_hz * self.PLANCK_EV_S
        self.e_c_ev = e_c_hz * self.PLANCK_EV_S
        
        # Delta in eV (convert from Hz if provided, else use material)
        if delta_l_hz is not None:
            self.delta_l_ev = delta_l_hz * self.PLANCK_EV_S
        else:
            self.delta_l_ev = self.delta_material
            
        if delta_r_hz is not None:
            self.delta_r_ev = delta_r_hz * self.PLANCK_EV_S
        else:
            self.delta_r_ev = self.delta_material
        
        # Computed quantities
        self.ej_ec_ratio = self.e_j_hz / self.e_c_hz
        self.curly_n = (self.dos * 
                       np.sqrt(2 * np.pi * self.delta_l_ev * 
                               self.KB_EV_K * self.temperature_k))
        
    @classmethod
    def from_capacitance(cls, c_total_f, delta_hz, r_n_ohm, **kwargs):
        """
        Create OCS instance from capacitance and resistance
        
        Parameters
        ----------
        c_total_f : float
            Total capacitance Cᵨ = Cⱼ + Cᵍ + Cₛₕᵤₙₜ [F]
        delta_hz : float
            Superconducting gap [Hz]
        r_n_ohm : float
            Normal state resistance [Ω]
        **kwargs
            Additional arguments passed to __init__
            
        Returns
        -------
        OCS
            Configured OCS instance
        """
        # Compute energies in eV first
        e_c_ev = cls.ELECTRON_CHARGE / (2 * c_total_f)  # eV
        delta_ev = delta_hz * cls.PLANCK_EV_S  # eV
        e_j_ev = (cls.PLANCK_EV_S * delta_ev / 
                 (8 * cls.ELECTRON_CHARGE * r_n_ohm))  # eV
        
        # Convert to Hz for the constructor
        e_c_hz = e_c_ev / cls.PLANCK_EV_S
        e_j_hz = e_j_ev / cls.PLANCK_EV_S
        
        return cls(e_j_hz, e_c_hz, r_n_ohm=r_n_ohm, **kwargs)
    
    def build_hamiltonian(self, offset_charge, charge_cutoff=18):
        """
        Construct Cooper-pair box Hamiltonian in charge basis
        
        The Hamiltonian is H = 4Eᴄ(n̂ - nᵍ)² - (Eⱼ/2)(|n⟩⟨n+1| + h.c.)
        where n̂ is the Cooper pair number operator and nᵍ is the 
        offset charge.
        
        Parameters
        ----------
        offset_charge : float
            Dimensionless offset charge nᵍ = CᵍVᵍ/(2e)
        charge_cutoff : int, optional
            Charge basis cutoff, default 18
            
        Returns
        -------
        h : ndarray
            Hamiltonian matrix [eV]
        """
        n_dim = 2 * charge_cutoff + 1
        h = np.zeros((n_dim, n_dim))
        
        # Charge states from -charge_cutoff to charge_cutoff
        charge_states = np.arange(-charge_cutoff, charge_cutoff + 1)
        
        # Diagonal: charging energy term 4Eᴄ(n - nᵍ)²
        h[np.arange(n_dim), np.arange(n_dim)] = (
            4.0 * self.e_c_ev * (charge_states - offset_charge) ** 2
        )
        
        # Off-diagonal: Josephson tunneling -Eⱼ/2
        h[np.arange(n_dim - 1), np.arange(1, n_dim)] = -self.e_j_ev / 2
        h[np.arange(1, n_dim), np.arange(n_dim - 1)] = -self.e_j_ev / 2
        
        return h
    
    def solve_eigensystem(self, offset_charge, charge_cutoff=18):
        """
        Diagonalize Hamiltonian to find energy eigenstates
        
        Parameters
        ----------
        offset_charge : float
            Dimensionless offset charge nᵍ
        charge_cutoff : int, optional
            Charge basis cutoff, default 18
            
        Returns
        -------
        eigenvalues : ndarray
            Energy eigenvalues [eV], sorted ascending
        eigenvectors : ndarray
            Energy eigenvectors as columns, sorted by eigenvalue
        """
        h = self.build_hamiltonian(offset_charge, charge_cutoff)
        eigenvalues, eigenvectors = eigh(h)
        return eigenvalues, eigenvectors
    
    def solve_system(self, offset_charges, num_levels=4, 
                    charge_cutoff=18):
        """
        Solve for even and odd parity energy levels
        
        Even parity: extra Cooper pair on island (nᵍ)
        Odd parity: extra unpaired electron (nᵍ + 0.5)
        
        Parameters
        ----------
        offset_charges : array_like
            Array of offset charge values to compute
        num_levels : int, optional
            Number of energy levels to return, default 4
        charge_cutoff : int, optional
            Charge basis cutoff, default 18
            
        Returns
        -------
        energies_even : ndarray
            Even parity energies [eV], shape (len(n_g), num_levels)
        energies_odd : ndarray
            Odd parity energies [eV], shape (len(n_g), num_levels)
        energy_diff : ndarray
            Parity splitting energy [eV], shape (len(n_g),)
        """
        offset_charges = np.atleast_1d(offset_charges)
        num_points = len(offset_charges)
        
        energies_even = np.zeros((num_points, num_levels))
        energies_odd = np.zeros((num_points, num_levels))
        
        for i, n_g in enumerate(offset_charges):
            # Even parity: integer charge
            evals_even, _ = self.solve_eigensystem(n_g, 
                                                   charge_cutoff)
            energies_even[i, :] = evals_even[:num_levels]
            
            # Odd parity: half-integer charge
            evals_odd, _ = self.solve_eigensystem(n_g + 0.5, 
                                                  charge_cutoff)
            energies_odd[i, :] = evals_odd[:num_levels]
        
        # Parity splitting includes superconducting gap difference
        energy_diff = (energies_odd[:, 0] - energies_even[:, 0] + 
                      self.delta_l_ev - self.delta_r_ev)
        
        return energies_even, energies_odd, energy_diff
    
    def compute_dispersive_matrix(self, offset_charge, coupling_g_hz, 
                                  readout_freq_hz, num_levels=6, parity='odd',
                                  charge_cutoff=30):
        """
        Compute dispersive shift and matrix elements
        
        The dispersive shift is:
        χᵢ = g² ∑ⱼ≠ᵢ [2ωᵢⱼ |⟨j,p|n̂|i,p⟩|² / (ωᵢⱼ² - ωᵣ²)]
        
        Parameters
        ----------
        offset_charge : float
            Dimensionless offset charge nᵍ
        coupling_g_hz : float
            Transmon-resonator coupling strength [Hz]
        readout_freq_hz : float
            Resonator frequency [Hz]
        num_levels : int, optional
            Number of levels for matrix elements, default 6
        parity : str, optional
            Parity of the qubit, default 'odd'. Options: 'odd', 'even'
        charge_cutoff : int, optional
            Charge basis cutoff (higher for accuracy), default 30
            
        Returns
        -------
        matrix_elements : ndarray
            Matrix elements |⟨j|n̂|0⟩|², shape (num_levels, num_levels)
        chi_ip : ndarray
            Dispersive shift χᵢ,ₚ [Hz], shape (num_levels,)
        """
        if parity == 'odd':
            offset_charge += 0.5
        elif parity == 'even':
            pass
        else:
            raise ValueError(f"Invalid parity: {parity}")
        
        # Solve eigensystem
        eigenvalues, eigenvectors = self.solve_eigensystem(
            offset_charge, charge_cutoff
        )
        
        # Build number operator n̂ = ∑ₙ (n - nᵍ)|n⟩⟨n|
        charge_states = np.arange(-charge_cutoff, charge_cutoff + 1)
        number_operator = np.diag(charge_states)
        
        # Transform to energy eigenbasis
        total_states = len(eigenvalues)
        matrix_elements = np.zeros((num_levels, num_levels))
        chi_ip = np.zeros(num_levels)
        
        for i in range(num_levels):
            for j in range(total_states):
                if i != j:
                    # Transition frequency
                    omega_ij = (eigenvalues[i] - eigenvalues[j]) / (
                        self.PLANCK_EV_S
                    )  # Hz
                    
                    # Matrix element |⟨j|n̂|i⟩|²
                    mat_elem_sq = np.abs(
                        eigenvectors[:, j].conj() @ number_operator @ 
                        eigenvectors[:, i]
                    ) ** 2
                    
                    # Store matrix elements for first num_levels states
                    if j < num_levels:
                        matrix_elements[i, j] = mat_elem_sq
                    
                    # Dispersive shift contribution
                    chi_contrib = (2.0 * omega_ij * mat_elem_sq / 
                                  (omega_ij ** 2 - readout_freq_hz ** 2))
                    chi_ip[i] += chi_contrib
        
        # Scale by coupling strength squared
        chi_ip *= coupling_g_hz ** 2
        
        return matrix_elements, chi_ip
    
    def plot_energy_levels(self, offset_charges=None, num_levels=5, 
                          freq_resonator_hz=None, coupling_g_hz=None,
                          figsize=(4, 3), ylim=None):
        """
        Plot energy level diagram vs offset charge
        
        Shows both even parity (solid) and odd parity (dashed) levels.
        
        Parameters
        ----------
        offset_charges : array_like, optional
            Offset charge values, default linspace(0, 1, 500)
        num_levels : int, optional
            Number of levels to plot, default 4
        figsize : tuple, optional
            Figure size (width, height), default (4, 3)
        ylim : tuple, optional
            Y-axis limits, default None
            
        Returns
        -------
        fig : matplotlib.figure.Figure
            Figure handle
        ax : matplotlib.axes.Axes
            Axes handle
        """
        if offset_charges is None:
            offset_charges = np.linspace(0, 1, 500)
        
        energies_even, energies_odd, energy_diff = (
            self.solve_system(offset_charges, num_levels)
        )
        
        # Convert to frequencies (GHz) relative to ground state
        freq_even = ((energies_even - energies_even[:, [0]]) / 
                    self.PLANCK_EV_S / 1e9) # GHz
        freq_odd = ((energies_odd - energies_odd[:, [0]]) / 
                   self.PLANCK_EV_S / 1e9) # GHz
        
        # Print key frequencies
        print("At offset charge 0:")
        print(f"Eⱼ/Eᴄ = {self.ej_ec_ratio:.2f}")
        print(f"E₀₁ = {freq_odd[0, 1]:.3f} GHz")
        print(f"E₀₂ = {freq_odd[0, 2]:.3f} GHz")
        print(f"E₀₃ = {freq_odd[0, 3]:.3f} GHz")
        print(f"ΔE/Eᴄ = {energy_diff[0] / self.e_c_ev:.4f}")
        
        with plt.style.context(self._style_path):
            fig, ax = plt.subplots(figsize=figsize)
            
            # Get colormaps
            cmap_even = cm.get_cmap('Reds')
            cmap_odd = cm.get_cmap('Blues')
            
            # Plot even parity (solid lines) - Reds colormap
            for j in range(num_levels):
                color = cmap_even(0.2 + 0.7 * j / max(num_levels - 1, 1))
                ax.plot(offset_charges, freq_even[:, j], linewidth=2, 
                       color=color,
                       label=f'|{j},e⟩')
            
            # Plot odd parity (dashed lines) - Blues colormap
            for j in range(num_levels):
                color = cmap_odd(0.3 + 0.7 * j / max(num_levels - 1, 1))
                ax.plot(offset_charges, freq_odd[:, j], 
                       linewidth=2, color=color,
                       label=f'|{j},o⟩')
            
            if freq_resonator_hz is not None:
                ax.axhline(freq_resonator_hz / 1e9, 
                          color='black', linestyle=':')

            ax.set_xlim([0, 1])
            ax.set_xlabel(r'Offset Charge [$C_g V_g / 2e$]')
            ax.set_ylabel(r'$f_{0i}$ [GHz]')
            
            # Construct comprehensive title
            title_parts = [
                f'$\\xi={self.ej_ec_ratio:.1f}$',
                f'$E_J={self.e_j_hz/1e9:.2f}$ GHz',
                f'$E_C={self.e_c_hz/1e9:.3f}$ GHz'
            ]
            if coupling_g_hz is not None:
                title_parts.append(f'$g={coupling_g_hz/1e6:.0f}$ MHz')
            #title_parts.append(f'$T={self.temperature_k*1e3:.0f}$ mK')
            ax.set_title(', '.join(title_parts), fontsize=7)
            
            ax.minorticks_on()
            if num_levels <= 4:
                ax.legend(loc="best")
            ax.grid(alpha=0.3)
            if ylim is not None:
                ax.set_ylim(ylim)
        
        return fig, ax
    
    def plot_matrix_elements(self, offset_charges=None, 
                            coupling_g_hz=150e6, 
                            readout_freq_hz=7.0e9, num_levels=6, parity='odd',
                            energy_level=0,
                            figsize=(4, 3)):
        """
        Plot charge matrix elements vs offset charge
        
        Shows |⟨j,o|n̂|0,o⟩| for transitions from ground state.
        
        Parameters
        ----------
        offset_charges : array_like, optional
            Offset charge values, default linspace(0, 1, 500)
        coupling_g_hz : float, optional
            Coupling strength [Hz], default 150 MHz
        readout_freq_hz : float, optional
            Resonator frequency [Hz], default 7 GHz
        num_levels : int, optional
            Number of levels, default 6
        parity : str, optional
            Parity of the qubit, default 'odd'. Options: 'odd', 'even'
        energy_level : int, optional
            Energy level to plot, default 0
        figsize : tuple, optional
            Figure size, default (4, 3)
            
        Returns
        -------
        fig, ax : matplotlib figure and axes
        """
        if offset_charges is None:
            offset_charges = np.linspace(0, 1, 500)
        
        num_points = len(offset_charges)
        matrix_elems = np.zeros((num_points, num_levels))
        
        for i, n_g in enumerate(offset_charges):
            mat, _ = self.compute_dispersive_matrix(
                n_g, coupling_g_hz, readout_freq_hz, num_levels, parity=parity
            )
            matrix_elems[i, :] = mat[energy_level, :]  # From ground state
        
        with plt.style.context(self._style_path):
            fig, ax = plt.subplots(figsize=figsize)
            
            # Get colormap for odd parity
            if parity == 'odd':
                cmap = cm.get_cmap('Blues')
            elif parity == 'even':
                cmap = cm.get_cmap('Reds')
            else:
                raise ValueError(f"Invalid parity: {parity}")
            
            # Plot only transitions (j>0)
            for j in range(1, num_levels):
                color = cmap(0.3 + 0.7 * j / max(num_levels - 1, 1))
                ax.semilogy(offset_charges, matrix_elems[:, j], 
                           linewidth=2, color=color, label=f'j={j}')
            
            ax.set_ylim([1e-5, 2e0])
            ax.set_xlim([0, 1])
            ax.set_xlabel(r'Offset Charge [$C_g V_g / 2e$]')
            ax.set_ylabel(rf'$|\langle j,{parity[0]}|\hat{{n}}|{energy_level},{parity[0]}\rangle|^2$')
            
            # Construct comprehensive title
            title_parts = [
                f'$\\xi={self.ej_ec_ratio:.1f}$',
                f'$E_J={self.e_j_hz/1e9:.2f}$ GHz',
                f'$E_C={self.e_c_hz/1e9:.3f}$ GHz',
                f'$g={coupling_g_hz/1e6:.0f}$ MHz',
                #f'$T={self.temperature_k*1e3:.0f}$ mK'
            ]
            ax.set_title(', '.join(title_parts), fontsize=7)
            
            ax.minorticks_on()
            ax.legend(loc="best")
            ax.grid(alpha=0.3, which='both')
        
        return fig, ax
    
    def plot_dispersive_shift(self, offset_charges=None, 
                             coupling_g_hz=150e6, 
                             readout_freq_hz=7.0e9,
                             num_levels=6, figsize=(4, 3),
                             ylim=[-10, 10]):
        """
        Plot dispersive shift χ vs offset charge for both parities
        
        Shows charge-parity-dependent shift for ground and first 
        excited states for both odd (blue) and even (red) parities.
        
        Parameters
        ----------
        offset_charges : array_like, optional
            Offset charge values, default linspace(0, 1, 500)
        coupling_g_hz : float, optional
            Coupling strength [Hz], default 150 MHz
        readout_freq_hz : float, optional
            Resonator frequency [Hz], default 7 GHz
        num_levels : int, optional
            Number of levels, default 6
        figsize : tuple, optional
            Figure size, default (4, 3)
        ylim : tuple, optional
            Y-axis limits, default [-10, 10] MHz
        Returns
        -------
        fig, ax : matplotlib figure and axes
        """
        if offset_charges is None:
            offset_charges = np.linspace(0, 1, 500)
        
        num_points = len(offset_charges)
        chi_vals_odd = np.zeros((num_points, num_levels))
        chi_vals_even = np.zeros((num_points, num_levels))
        
        # Compute chi for both parities
        for i, n_g in enumerate(offset_charges):
            _, chi_odd = self.compute_dispersive_matrix(
                n_g, coupling_g_hz, readout_freq_hz, num_levels, 
                parity='odd'
            )
            _, chi_even = self.compute_dispersive_matrix(
                n_g, coupling_g_hz, readout_freq_hz, num_levels, 
                parity='even'
            )
            chi_vals_odd[i, :] = chi_odd
            chi_vals_even[i, :] = chi_even
        
        # Compute f_10 at sweet spot (n_g = 0.5) for reference lines
        eigenvalues, _ = self.solve_eigensystem(0.5, charge_cutoff=30)
        f_10 = (eigenvalues[1] - eigenvalues[0]) / self.PLANCK_EV_S
        
        # Compute approximate chi formulas
        chi_approx_1 = - (coupling_g_hz ** 2) / (f_10-readout_freq_hz)  # Hz
        chi_approx_2 = - (coupling_g_hz ** 2) / (f_10-readout_freq_hz) + \
                       (coupling_g_hz ** 2) / ((f_10-readout_freq_hz) - self.e_c_hz)  # Hz
        
        with plt.style.context(self._style_path):
            fig, ax = plt.subplots(figsize=figsize)
            
            # Get colormaps for both parities
            cmap_odd = cm.get_cmap('Blues')
            cmap_even = cm.get_cmap('Reds')
            
            # Plot both parities
            for j in range(num_levels-2):
                # Odd parity (blue)
                color_odd = cmap_odd(0.3 + 0.7 * j / max(num_levels - 3, 1))
                ax.plot(offset_charges, chi_vals_odd[:, j] / 1e6, 
                       linewidth=2, color=color_odd, 
                       label=f'|{j},o⟩')
                
                # Even parity (red)
                color_even = cmap_even(0.3 + 0.7 * j / max(num_levels-3, 1))
                ax.plot(offset_charges, chi_vals_even[:, j] / 1e6, 
                       linewidth=2, color=color_even, 
                       label=f'|{j},e⟩')
            
            ax.set_xlim([0, 1])
            ax.set_ylim(ylim)
            ax.set_xlabel(r'Offset Charge [$C_g V_g / 2e$]')
            ax.set_ylabel(r'$\chi_{i,p}$ [MHz]')
            
            # Add approximate formula reference lines
            ax.axhline(chi_approx_1 / 1e6, color='gray', 
                      linestyle='--', linewidth=1.5, alpha=0.7,
                      label=r'$-g^2/\Delta$')
            ax.axhline(chi_approx_2 / 1e6, color='black', 
                      linestyle=':', linewidth=1.5, alpha=0.7,
                      label=r'$-g^2/\Delta+g^2/(\Delta-E_C)$')
            
            # Construct comprehensive title
            title_parts = [
                f'$\\xi={self.ej_ec_ratio:.1f}$',
                f'$E_J={self.e_j_hz/1e9:.2f}$ GHz',
                f'$E_C={self.e_c_hz/1e9:.3f}$ GHz',
                f'$g={coupling_g_hz/1e6:.0f}$ MHz',
                f'$f_r={readout_freq_hz/1e9:.2f}$ GHz',
                #f'$T={self.temperature_k*1e3:.0f}$ mK'
            ]
            ax.set_title(', '.join(title_parts), fontsize=7)
            
            ax.minorticks_on()
            ax.legend(loc="best", ncol=2, fontsize=6)
            ax.grid(alpha=0.3)
        
        return fig, ax
    
    def plot_parity_shift_vs_frequency(self, freq_range_hz=None, 
                                      coupling_g_hz=150e6, 
                                      num_levels=6,
                                      freq_min_hz=1.0e9,
                                      offset_charges=[0.5],
                                      figsize=(4, 3)):
        """
        Plot parity-dependent dispersive shift vs resonator frequency
        
        Plots the difference between odd and even parity chi values
        for multiple offset charges.
        
        Parameters
        ----------
        freq_range_hz : array_like, optional
            Resonator frequencies [Hz], auto-computed if None
        coupling_g_hz : float, optional
            Coupling strength [Hz], default 150 MHz
        num_levels : int, optional
            Number of levels, default 6
        freq_min_hz : float, optional
            Minimum frequency for range, default 1 GHz
        offset_charges : array_like, optional
            Offset charge values to plot, default [0.5]
        figsize : tuple, optional
            Figure size, default (4, 3)
            
        Returns
        -------
        fig : matplotlib figure
        (ax1, ax2) : tuple of matplotlib axes
            ax1: top subplot with parity shift
            ax2: bottom subplot with detuning ratio
        """
        # Auto-determine frequency range if not provided
        if freq_range_hz is None:
            _, energies_odd, _ = self.solve_system([0.5], 4)
            freq_min = freq_min_hz
            freq_odd_3 = ((energies_odd[0, 3] - energies_odd[0, 0]) / 
                         self.PLANCK_EV_S)  # f_03
            freq_max = freq_odd_3 + 0.2e9  # Add 200 MHz
            freq_range_hz = np.linspace(freq_min, freq_max, 500)
        
        # Calculate transition frequencies for markers
        _, energies, _ = self.solve_system([0], num_levels)
        freq_10_0 = ((energies[0, 1] - energies[0, 0]) / 
                  self.PLANCK_EV_S)
        freq_20_0 = ((energies[0, 2] - energies[0, 0]) / 
                  self.PLANCK_EV_S)
        freq_30_0 = ((energies[0, 3] - energies[0, 0]) / 
                  self.PLANCK_EV_S)
        _, energies, _ = self.solve_system([0.5], num_levels)
        freq_10_05 = ((energies[0, 1] - energies[0, 0]) / 
                  self.PLANCK_EV_S)
        freq_20_05 = ((energies[0, 2] - energies[0, 0]) / 
                  self.PLANCK_EV_S)
        freq_30_05 = ((energies[0, 3] - energies[0, 0]) / 
                  self.PLANCK_EV_S)
        # Median frequency
        freq_10 = (freq_10_0 + freq_10_05) / 2
        freq_20 = (freq_20_0 + freq_20_05) / 2
        freq_30 = (freq_30_0 + freq_30_05) / 2

        # Compute chi_diff for each offset charge
        # Result shape: (len(freq_range), len(offset_charges))
        # Transpose to get: (len(offset_charges), len(freq_range))
        chi_diffs = self.compute_delta_chi_0(
            freq_range_hz, offset_charges, coupling_g_hz, num_levels
        ).T
        
        with plt.style.context(self._style_path):
            # Create 2 subplots with shared x-axis
            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=figsize, 
                                          sharex=True,
                                          gridspec_kw={'height_ratios': [3, 1],
                                                      'hspace': 0.01})
            
            # === Top subplot: parity shift ===
            # Mark f10 region with axvspan
            ax1.axvspan(freq_10_0 / 1e9, freq_10_05 / 1e9,
                       alpha=0.2, color='gray', label=r'$f_{10}$', lw=0)
            
            # Mark f10 ± g region with axvspan
            ax1.axvspan((freq_10 - coupling_g_hz) / 1e9, 
                       (freq_10 + coupling_g_hz) / 1e9,
                       alpha=0.15, color='red', 
                       label=r'$f_{10} \pm g$', lw=0)
            
            # Mark f20 and f30 regions with axvspan
            ax1.axvspan(freq_20_0 / 1e9, freq_20_05 / 1e9,
                       alpha=0.2, color='C0', label=r'$f_{20}$', lw=0)
            ax1.axvspan(freq_30_0 / 1e9, freq_30_05 / 1e9,
                       alpha=0.2, color='C1', label=r'$f_{30}$', lw=0)
            
            # Get colormap
            cmap = cm.get_cmap('magma')
            
            # Plot curves for each offset charge
            for idx, n_g in enumerate(offset_charges):
                color = cmap(idx / max(len(offset_charges) - 1, 1))
                ax1.semilogy(freq_range_hz / 1e9, 
                           np.abs(chi_diffs[idx, :]) / 1e6, 
                           linewidth=2, color=color, 
                           label=f'$n_g={n_g:.2f}$')
            
            ax1.set_ylabel(r'$|\chi_{0, o}-\chi_{0, e}|$ [MHz]')
            
            # Construct comprehensive title
            title_parts = [
                f'$\\xi={self.ej_ec_ratio:.1f}$',
                f'$E_J={self.e_j_hz/1e9:.2f}$ GHz',
                f'$E_C={self.e_c_hz/1e9:.3f}$ GHz',
                f'$g={coupling_g_hz/1e6:.0f}$ MHz',
                #f'$T={self.temperature_k*1e3:.0f}$ mK'
            ]
            ax1.set_title(', '.join(title_parts), fontsize=7)
            ax1.set_xlim(freq_range_hz[0] / 1e9, freq_range_hz[-1] / 1e9)
            ax1.set_ylim(1e-3)
            ax1.minorticks_on()
            ax1.legend(loc="best", fontsize=6, ncol=3)
            ax1.grid(alpha=0.3, which='both')
            
            # === Bottom subplot: detuning ratio ===
            # Mark f10 region with axvspan
            ax2.axvspan(freq_10_0 / 1e9, freq_10_05 / 1e9,
                       alpha=0.2, color='gray', lw=0)
            
            # Mark f10 ± g region with axvspan
            ax2.axvspan((freq_10 - coupling_g_hz) / 1e9, 
                       (freq_10 + coupling_g_hz) / 1e9,
                       alpha=0.15, color='red', lw=0)
            
            # Mark f20 and f30 regions with axvspan
            ax2.axvspan(freq_20_0 / 1e9, freq_20_05 / 1e9,
                       alpha=0.2, color='C0', lw=0)
            ax2.axvspan(freq_30_0 / 1e9, freq_30_05 / 1e9,
                       alpha=0.2, color='C1', lw=0)
            
            # Calculate detuning ratio and mask region between f_10_0 and f_10_05
            detuning_ratio = (np.abs(freq_range_hz - freq_10) / 
                            np.abs(freq_range_hz + freq_10))
            # Set values between f_10_0 and f_10_05 to NaN to create a gap
            mask_f10 = ((freq_range_hz >= min(freq_10_0, freq_10_05)) & 
                       (freq_range_hz <= max(freq_10_0, freq_10_05)))
            detuning_ratio[mask_f10] = np.nan
            
            ax2.semilogy(freq_range_hz / 1e9, detuning_ratio, 
                        linewidth=1, color='black')
            
            ax2.set_xlim(freq_range_hz[0] / 1e9, freq_range_hz[-1] / 1e9)
            ax2.set_xlabel('Readout Frequency [GHz]')
            ax2.set_ylabel(r'$|f_r - f_{10}|/|f_r+f_{10}|$')
            ax2.minorticks_on()
            ax2.grid(alpha=0.3, which='both')
        
        # Print summary
        freq_02 = ((energies[0, 2] - energies[0, 0]) / 
                  self.PLANCK_EV_S)
        anharmonicity = freq_02 - 2 * freq_10
        
        # Chi at resonator for qubit state readout (simple estimate)
        if hasattr(self, '_readout_freq_hz'):
            f_r_use = self._readout_freq_hz
        else:
            f_r_use = 7.0e9
        chi_resonator = (coupling_g_hz ** 2 * anharmonicity / 
                        np.abs(f_r_use - freq_10) / 
                        (np.abs(f_r_use - freq_10) + anharmonicity))
        
        print(f"\nf_10: {freq_10 / 1e9:.3f} GHz")
        print(f"f_20: {freq_20 / 1e9:.3f} GHz")
        print(f"f_30: {freq_30 / 1e9:.3f} GHz")
        print(f"Resonator frequency: {f_r_use / 1e9:.2f} GHz")
        print(f"χ (resonator state shift): {chi_resonator / 1e6:.3f} MHz")
        
        return fig, (ax1, ax2)
    
    def plot_parity_shift_vs_ng(self, offset_charges=None, 
                                coupling_g_hz=150e6, 
                                num_levels=6,
                                readout_freqs=[7.0e9],
                                figsize=(4, 3)):
        """
        Plot parity-dependent dispersive shift vs offset charge
        
        Plots the difference between odd and even parity chi values
        as a function of offset charge for multiple resonator frequencies.
        
        Parameters
        ----------
        offset_charges : array_like, optional
            Offset charge values from 0 to 1, default linspace(0, 1, 500)
        coupling_g_hz : float, optional
            Coupling strength [Hz], default 150 MHz
        num_levels : int, optional
            Number of levels, default 6
        readout_freqs : array_like, optional
            Resonator frequencies [Hz] to plot, default [7.0e9]
        figsize : tuple, optional
            Figure size, default (4, 3)
            
        Returns
        -------
        fig, ax : matplotlib figure and axes
        """
        if offset_charges is None:
            offset_charges = np.linspace(0, 1, 500)
        
        # Calculate f10 for reference
        _, energies, _ = self.solve_system([0], num_levels)
        freq_10 = ((energies[0, 1] - energies[0, 0]) / 
                  self.PLANCK_EV_S)
        
        # Compute chi_diff for each resonator frequency
        chi_diffs = self.compute_delta_chi_0(
            readout_freqs, offset_charges, coupling_g_hz, num_levels
        )
        
        with plt.style.context(self._style_path):
            fig, ax = plt.subplots(figsize=figsize)
            
            # Get colormap
            cmap = cm.get_cmap('magma')
            
            # Plot curves for each resonator frequency
            for idx, f_r in enumerate(readout_freqs):
                color = cmap(idx / max(len(readout_freqs) - 1, 1))
                ax.semilogy(offset_charges, 
                           np.abs(chi_diffs[idx, :]) / 1e6, 
                           linewidth=2, color=color, 
                           label=f'$f_r={f_r/1e9:.1f}$ GHz')
            
            ax.set_xlim([0, 1])
            ax.set_ylim(1e-3)
            ax.set_xlabel(r'Offset Charge [$C_g V_g / 2e$]')
            ax.set_ylabel(r'$|\chi_{0, o} - \chi_{0, e}|$ [MHz]')
            
            # Construct comprehensive title
            title_parts = [
                f'$\\xi={self.ej_ec_ratio:.1f}$',
                f'$E_J={self.e_j_hz/1e9:.2f}$ GHz',
                f'$E_C={self.e_c_hz/1e9:.3f}$ GHz',
                f'$g={coupling_g_hz/1e6:.0f}$ MHz',
                #f'$T={self.temperature_k*1e3:.0f}$ mK'
            ]
            ax.set_title(', '.join(title_parts), fontsize=7)
            
            ax.minorticks_on()
            ax.legend(loc="best", fontsize=7)
            ax.grid(alpha=0.3, which='both')
        
        # Print summary
        print(f"\nf_10: {freq_10 / 1e9:.3f} GHz")
        
        return fig, ax
    
    def compute_average_chi_0(self, coupling_g_hz, readout_freq_hz, 
                             num_points=500, num_levels=6):
        """
        Compute average dispersive shift for ground state (i=0)
        
        Averages chi_0 over offset charges from 0 to 1 for both parities.
        
        Parameters
        ----------
        coupling_g_hz : float
            Coupling strength [Hz]
        readout_freq_hz : float
            Resonator frequency [Hz]
        num_points : int, optional
            Number of offset charge points to sample, default 500
        num_levels : int, optional
            Number of levels for dispersive calculation, default 6
            
        Returns
        -------
        chi_0_avg : float
            Average dispersive shift for ground state [Hz]
        """
        offset_charges = np.linspace(0, 1, num_points)
        chi_0_odd = np.zeros(num_points)
        chi_0_even = np.zeros(num_points)
        
        for i, n_g in enumerate(offset_charges):
            # Compute for odd parity
            _, chi_odd = self.compute_dispersive_matrix(
                n_g, coupling_g_hz, readout_freq_hz, num_levels, 
                parity='odd'
            )
            chi_0_odd[i] = chi_odd[0]
            
            # Compute for even parity
            _, chi_even = self.compute_dispersive_matrix(
                n_g, coupling_g_hz, readout_freq_hz, num_levels, 
                parity='even'
            )
            chi_0_even[i] = chi_even[0]
        
        # Average over offset charges and both parities
        chi_0_odd_avg = np.mean(chi_0_odd)
        chi_0_even_avg = np.mean(chi_0_even)
        chi_0_avg = (chi_0_odd_avg + chi_0_even_avg) / 2
        
        return chi_0_avg
    
    def compute_delta_chi_0(self, readout_freqs, offset_charges, 
                           coupling_g_hz, num_levels=6):
        """
        Compute parity-dependent dispersive shift difference
        
        Computes Δχ₀ = χ₀,odd - χ₀,even for the ground state across
        different resonator frequencies and offset charges.
        
        Parameters
        ----------
        readout_freqs : array_like
            Resonator frequencies [Hz]
        offset_charges : array_like
            Offset charge values (dimensionless)
        coupling_g_hz : float
            Coupling strength [Hz]
        num_levels : int, optional
            Number of levels for dispersive calculation, default 6
            
        Returns
        -------
        chi_diffs : ndarray
            Parity shift difference [Hz], 
            shape (len(readout_freqs), len(offset_charges))
        """
        readout_freqs = np.atleast_1d(readout_freqs)
        offset_charges = np.atleast_1d(offset_charges)
        
        # Store chi_diff for each resonator frequency
        chi_diffs = np.zeros((len(readout_freqs), len(offset_charges)))
        
        for idx, f_r in enumerate(readout_freqs):
            for i, n_g in enumerate(offset_charges):
                # Chi for odd parity
                _, chi_odd = self.compute_dispersive_matrix(
                    n_g, coupling_g_hz, f_r, num_levels, 
                    parity='odd'
                )
                # Chi for even parity
                _, chi_even = self.compute_dispersive_matrix(
                    n_g, coupling_g_hz, f_r, num_levels, 
                    parity='even'
                )
                chi_diffs[idx, i] = chi_odd[0] - chi_even[0]
        
        return chi_diffs
    
    def estimate_g_from_chi_0(self, chi_0_measured_hz, readout_freq_hz, 
                             g_initial_hz=100e6, num_levels=6):
        """
        Estimate coupling g from target average dispersive shift chi_0
        
        Uses iterative root-finding to determine g that produces the 
        desired average dispersive shift.
        
        Parameters
        ----------
        chi_0_measured_hz : float
            Measured dispersive shift by EPR Chi Matrix (ND) [Hz]
        readout_freq_hz : float
            Resonator frequency [Hz]
        g_initial_hz : float, optional
            Initial guess for coupling [Hz], default 100 MHz
        num_levels : int, optional
            Number of levels for calculation, default 6
            
        Returns
        -------
        g_estimated_hz : float
            Estimated coupling strength [Hz]
        chi_0_final_hz : float
            Final achieved average dispersive shift [Hz]
        """
        from scipy.optimize import fsolve
        
        def objective(g):
            """Objective function: chi_0(g) - target"""
            chi_0_avg = self.compute_average_chi_0(
                g, readout_freq_hz, num_points=100, 
                num_levels=num_levels
            )
            return chi_0_avg - chi_0_measured_hz
        
        # Solve for g
        g_estimated_hz = fsolve(objective, g_initial_hz)[0]
        
        return g_estimated_hz
    
    def plot_all(self, offset_charges=None, coupling_g_hz=150e6, 
                readout_freq_hz=7.0e9, num_levels=5):
        """
        Generate all standard plots for OCS transmon analysis
        
        Creates four plots:
        1. Energy level diagram
        2. Matrix elements
        3. Dispersive shift vs offset charge
        4. Parity shift vs resonator frequency
        
        Parameters
        ----------
        offset_charges : array_like, optional
            Offset charge values
        coupling_g_hz : float, optional
            Coupling strength [Hz], default 150 MHz
        readout_freq_hz : float, optional
            Resonator frequency [Hz], default 7 GHz
        num_levels : int, optional
            Number of levels, default 6
            
        Returns
        -------
        figs : list
            List of figure handles
        """
        self._readout_freq_hz = readout_freq_hz  # Store for later use
        
        print("=" * 60)
        print(f"OCS Transmon Parameters:")
        print(f"Material: {self.material_name} (Tc = {self.tc:.3f} K)")
        print(f"Eⱼ = {self.e_j_ev / self.KB_EV_K:.3f} K·kᴮ = "
              f"{self.e_j_hz / 1e9:.3f} GHz")
        print(f"Eᴄ = {self.e_c_ev / self.KB_EV_K:.4f} K·kᴮ = "
              f"{self.e_c_hz / 1e9:.4f} GHz")
        print(f"Eⱼ/Eᴄ = {self.ej_ec_ratio:.2f}")
        print(f"Rₙ = {self.r_n_ohm / 1e3:.1f} kΩ")
        print(f"Δ = {self.delta_material / self.KB_EV_K * 1e3:.3f} mK·kᴮ")
        print(f"Resonator frequency = {readout_freq_hz / 1e9:.2f} GHz")
        print(f"Coupling g = {coupling_g_hz / 1e6:.1f} MHz")
        print(f"FWHM = {readout_freq_hz / 10000 / 1e6:.2f} MHz")
        print("=" * 60)
        
        figs = []
        
        # Figure 1: Energy levels
        print("\n[1/4] Plotting energy levels...")
        fig1, _ = self.plot_energy_levels(
            offset_charges, num_levels, readout_freq_hz, coupling_g_hz
        )
        figs.append(fig1)
        fig1.show()
        
        # Figure 2: Matrix elements
        print("[2/4] Plotting matrix elements...")
        fig2, _ = self.plot_matrix_elements(
            None, coupling_g_hz, readout_freq_hz, num_levels
        )
        figs.append(fig2)
        fig2.show()
            
        # Figure 3: Dispersive shift
        print("[3/4] Plotting dispersive shift...")
        fig3, _ = self.plot_dispersive_shift(
            None, coupling_g_hz, readout_freq_hz, num_levels, ylim=[-20, 20]
        )
        figs.append(fig3)
        fig3.show()
        
        # Figure 4: Parity shift vs frequency
        print("[4/4] Plotting parity shift vs frequency...")
        fig4, _ = self.plot_parity_shift_vs_frequency(
            None, coupling_g_hz, num_levels
        )
        figs.append(fig4)
        fig4.show()
        
        print("\nAll plots complete!")
        return figs


def main():
    """Example usage of OCS class"""
    # Example from the MATLAB script (WashU parameters)
    ej_ec_ratio = 12
    e_j_hz = 8.335e9  # ~8.3 GHz (0.4 K·kB)
    e_c_hz = e_j_hz / ej_ec_ratio  # Hz
    
    # Create OCS instance
    ocs = OCS(
        e_j_hz=e_j_hz,
        e_c_hz=e_c_hz,
        temperature_k=0.012,  # 12 mK
        r_n_ohm=27e3  # 27 kΩ
    )
    
    # Generate all plots
    figs = ocs.plot_all(
        coupling_g_hz=150e6,  # 150 MHz
        readout_freq_hz=7.0e9,  # 7 GHz
        num_levels=5
    )
    
    plt.show()
    
    return ocs, figs


if __name__ == "__main__":
    ocs, figs = main()

