"""
QPD Transmon Analysis

This module implements simulation and analysis tools for QPD transmons,
based on the physics described in:
- Serniak et al., PRA (2019): https://arxiv.org/pdf/1903.00113
- Additional context from https://arxiv.org/pdf/2405.17192

The QPD transmon is in the intermediate regime between Cooper-pair box 
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


class QPD:
    """
    QPD Transmon Simulator
    
    This class simulates QPD transmons by diagonalizing the Cooper-pair
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
    _style_path = Path(__file__).parent / "qpd.mplstyle"
    
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
        Initialize QPD transmon simulator
        
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
        Create QPD instance from capacitance and resistance
        
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
        QPD
            Configured QPD instance
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

    def compute_quantum_capacitance(self, offset_charges,
                                    c_g_f=None, charge_cutoff=18,
                                    level=0):
        """
        Compute the quantum capacitance C_Q(n_g) for both parities.

        Numerically differentiates the parity-resolved level energy twice
        with respect to the dimensionless offset charge n_g (units of 2e):

            C_Q[F] = -(C_g / 2e)^2 · h · ∂²E_level[Hz] / ∂n_g²

        If `c_g_f` is None, the prefactor is dropped and the function
        instead returns the intrinsic second derivative ∂²E[Hz]/∂n_g²
        in Hz so callers can apply their own scaling (this is what the
        fitter consumes — the absolute scale is fit as a free parameter).

        Parameters
        ----------
        offset_charges : array_like
            Dimensionless offset-charge grid n_g where C_Q is evaluated.
            Should be evenly spaced; finite-difference accuracy degrades
            otherwise.
        c_g_f : float, optional
            Gate capacitance C_g [F]. If provided, the result is in Farads.
            If None (default), the result is in Hz (intrinsic ∂²E/∂n_g²/h).
        charge_cutoff : int, optional
            Charge basis cutoff passed through to `solve_system`,
            default 18.
        level : int, optional
            Energy level whose curvature is taken (0 = ground state).
            Default 0.

        Returns
        -------
        cq_even : ndarray
            C_Q for even parity at each n_g, shape (len(offset_charges),).
        cq_odd : ndarray
            C_Q for odd parity at each n_g, shape (len(offset_charges),).

        Notes
        -----
        - C_Q is computed via two successive `np.gradient` calls. To
          avoid the one-sided-stencil edge artifact at the boundaries
          of the requested grid, the inner solve is run on a grid
          padded by 2 points on each side; the padded points are then
          discarded so the returned arrays match the requested length.
          This requires a few extra eigensolves but keeps the central-
          difference accuracy uniform across the grid.
        - Sign: peaks of C_Q (charge dispersion) appear at n_g = 0.5 for
          even parity and n_g = 0 (or 1) for odd parity, consistent with
          the offset_charge → offset_charge + 0.5 convention used by
          `solve_system`.
        """
        offset_charges = np.atleast_1d(np.asarray(offset_charges,
                                                  dtype=float))

        if offset_charges.size < 2:
            raise ValueError(
                "compute_quantum_capacitance needs at least 2 grid points"
            )

        # Pad both ends of the grid so np.gradient can use central
        # differences everywhere we report.
        pad = 2
        dn = offset_charges[1] - offset_charges[0]
        left = offset_charges[0] - dn * np.arange(pad, 0, -1)
        right = offset_charges[-1] + dn * np.arange(1, pad + 1)
        ng_ext = np.concatenate([left, offset_charges, right])

        num_levels = max(level + 1, 2)
        e_even_ev, e_odd_ev, _ = self.solve_system(
            ng_ext, num_levels=num_levels,
            charge_cutoff=charge_cutoff,
        )

        # Convert from eV to Hz once, then differentiate on the
        # padded grid and trim back to the requested range.
        e_even_hz = e_even_ev[:, level] / self.PLANCK_EV_S
        e_odd_hz = e_odd_ev[:, level] / self.PLANCK_EV_S

        d2_even_hz = np.gradient(np.gradient(e_even_hz, ng_ext),
                                 ng_ext)[pad:-pad]
        d2_odd_hz = np.gradient(np.gradient(e_odd_hz, ng_ext),
                                ng_ext)[pad:-pad]

        if c_g_f is None:
            # Intrinsic curvature in Hz (sign convention matches C_Q in F).
            return -d2_even_hz, -d2_odd_hz

        prefactor = (c_g_f / (2.0 * self.ELECTRON_CHARGE)) ** 2 * h
        cq_even = -prefactor * d2_even_hz
        cq_odd = -prefactor * d2_odd_hz
        return cq_even, cq_odd

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

    def fit_quantum_capacitance(self, offset_charges, cq_measured,
                                parity='odd', p0=None,
                                charge_cutoff=18, fit_offset=True,
                                fit_scale=True, fixed_scale=1.0,
                                fit_baseline=False, fixed_baseline=0.0,
                                e_j_bounds_hz=(0.01e9, 100e9),
                                e_c_bounds_hz=(0.01e9, 100e9),
                                n_g0_bounds=(-0.5, 0.5),
                                scale_bounds=(-np.inf, np.inf),
                                baseline_bounds=(-np.inf, np.inf),
                                cq_sigma=None, print_level=0):
        """
        Fit a measured C_Q(n_g) trace to extract E_J, E_C, n_g0, scale,
        and optionally a constant y-offset `baseline`.

        Uses iminuit (MIGRAD + HESSE) on a least-squares cost. The model
        evaluates the *intrinsic* second derivative -∂²E[Hz]/∂n_g² for
        the requested parity using the numeric Hamiltonian (same
        primitive as `compute_quantum_capacitance`). The full model is

            cq_meas ≈ scale · f_intrinsic(n_g − n_g0; E_J, E_C) + baseline

        where the multiplicative `scale` absorbs the prefactor relating
        the dimensionless intrinsic curvature to whatever units the user
        supplied (Farads if cq_measured is in F, etc.), so a known gate
        capacitance is *not* required, and the additive `baseline`
        absorbs any constant DC offset in the readout chain.

        Parameters
        ----------
        offset_charges : array_like
            Dimensionless n_g grid where the trace was sampled.
        cq_measured : array_like
            Measured C_Q values at the corresponding n_g.
        parity : {'odd', 'even'}, optional
            Which parity branch to fit, default 'odd'.
        p0 : sequence, optional
            Initial guess (e_j_hz, e_c_hz, n_g0, scale). Defaults
            (self.e_j_hz, self.e_c_hz, 0.0, max|C_Q|/max|d²E/dn_g²|).
        charge_cutoff : int, optional
            Charge-basis cutoff for the inner Hamiltonian solve.
        fit_offset : bool, optional
            If False, fix n_g0 = 0 and only fit (E_J, E_C, scale).
        fit_scale : bool, optional
            If True (default) fit the overall amplitude `scale` jointly
            with the energies. If the data has been pre-converted to a
            known absolute convention (e.g. intrinsic ∂²E/∂n_g² in Hz),
            set `fit_scale=False` and pass `fixed_scale` — this breaks
            the (E_C, scale) shape/amplitude degeneracy and lets the
            ratio fit be tight.
        fixed_scale : float, optional
            Scale used when `fit_scale=False`. Default 1.0.
        fit_baseline : bool, optional
            If True, fit an additive constant `baseline` that absorbs
            an unknown DC offset on cq_measured. Default False (no
            offset, i.e. `baseline = fixed_baseline = 0`). The
            baseline is robustly identified when fit jointly with
            `scale`, but in the deep-transmon limit (E_J/E_C ≳ 10)
            f_intrinsic becomes nearly sinusoidal and the DC level is
            correlated with scale × charge-dispersion(E_J, E_C); the
            (E_J, E_C) ratio can then drift more than in a
            baseline-fixed fit. If you have an independent baseline
            estimate, prefer `fit_baseline=False, fixed_baseline=...`.
        fixed_baseline : float, optional
            Baseline used when `fit_baseline=False`. Default 0.0.
        e_j_bounds_hz, e_c_bounds_hz : tuple of (lo, hi), optional
            Hard bounds on E_J and E_C in Hz. Default (0.01 GHz,
            100 GHz) for both — covers essentially all useful
            experimental QPD/transmon parameter space.
        n_g0_bounds, scale_bounds, baseline_bounds : tuple of (lo, hi)
            Bounds on n_g0 (default ±0.5), scale (default unbounded),
            and baseline (default unbounded). Pass `(None, None)` to
            remove a bound.
        cq_sigma : float or array_like, optional
            Per-point uncertainty on cq_measured used to weight the
            least-squares cost. Defaults to a constant value derived
            from the data magnitude (does not affect the central
            values; only the reported errors and chi² scale).
        print_level : int, optional
            iminuit verbosity (0 silent, 1 short, 2 detailed). Default 0.

        Returns
        -------
        result : dict
            Keys: 'e_j_hz', 'e_c_hz', 'n_g0', 'scale', 'baseline',
            'ej_ec_ratio',
            'errors' (1σ HESSE uncertainties for each parameter,
            plus 'ej_ec_ratio' propagated from the (E_J, E_C)
            covariance),
            'cq_fit' (model evaluated at offset_charges),
            'residuals' (cq_measured - cq_fit),
            'fval' (final chi² value),
            'minuit' (the underlying iminuit.Minuit instance, for
            inspection of covariance, MINOS errors, contours, etc.).
        """
        from iminuit import Minuit

        offset_charges = np.asarray(offset_charges, dtype=float)
        cq_measured = np.asarray(cq_measured, dtype=float)

        if parity not in ('odd', 'even'):
            raise ValueError(f"Invalid parity: {parity!r}")

        if cq_sigma is None:
            sigma = max(np.std(cq_measured) * 1e-2,
                        np.max(np.abs(cq_measured)) * 1e-3)
            sigma = np.full_like(cq_measured, sigma)
        else:
            sigma = np.broadcast_to(np.asarray(cq_sigma, dtype=float),
                                    cq_measured.shape).copy()

        def _intrinsic_cq(ng_grid, e_j_hz, e_c_hz):
            tmp = QPD.__new__(QPD)
            tmp.e_j_hz = e_j_hz
            tmp.e_c_hz = e_c_hz
            tmp.e_j_ev = e_j_hz * QPD.PLANCK_EV_S
            tmp.e_c_ev = e_c_hz * QPD.PLANCK_EV_S
            tmp.delta_l_ev = 0.0
            tmp.delta_r_ev = 0.0
            cq_e, cq_o = QPD.compute_quantum_capacitance(
                tmp, ng_grid, c_g_f=None, charge_cutoff=charge_cutoff,
            )
            return cq_o if parity == 'odd' else cq_e

        def model_eval(e_j_hz, e_c_hz, n_g0, scale, baseline):
            return scale * _intrinsic_cq(offset_charges - n_g0,
                                         e_j_hz, e_c_hz) + baseline

        def cost(e_j_hz, e_c_hz, n_g0, scale, baseline):
            resid = (model_eval(e_j_hz, e_c_hz, n_g0, scale, baseline)
                     - cq_measured) / sigma
            return float(np.sum(resid * resid))

        if p0 is None:
            cq0 = _intrinsic_cq(offset_charges, self.e_j_hz, self.e_c_hz)
            denom = np.max(np.abs(cq0))
            if fit_baseline:
                # Use peak-to-peak amplitude (baseline absorbs the DC),
                # and seed baseline at <cq_measured> − scale·<cq0>.
                amp_data = 0.5 * (np.max(cq_measured)
                                  - np.min(cq_measured))
                scale_guess = (amp_data / denom if denom > 0 else 1.0)
                baseline_guess = float(np.mean(cq_measured)
                                       - scale_guess * np.mean(cq0))
            else:
                amp_data = np.max(np.abs(cq_measured))
                scale_guess = (amp_data / denom if denom > 0 else 1.0)
                baseline_guess = fixed_baseline
            ej0, ec0, ng0_0 = self.e_j_hz, self.e_c_hz, 0.0
        else:
            p0 = list(p0)
            ej0, ec0 = p0[0], p0[1]
            idx = 2
            ng0_0 = p0[idx] if (fit_offset and len(p0) > idx) else 0.0
            if fit_offset:
                idx += 1
            scale_guess = (p0[idx] if (fit_scale and len(p0) > idx)
                           else fixed_scale)
            if fit_scale:
                idx += 1
            baseline_guess = (p0[idx] if (fit_baseline and len(p0) > idx)
                              else fixed_baseline)

        scale0 = scale_guess if fit_scale else fixed_scale
        baseline0 = baseline_guess if fit_baseline else fixed_baseline

        m = Minuit(cost, e_j_hz=ej0, e_c_hz=ec0,
                   n_g0=ng0_0, scale=scale0, baseline=baseline0)
        m.errordef = Minuit.LEAST_SQUARES
        m.print_level = print_level

        # Step sizes: O(parameter magnitude) for energies, O(0.01) for
        # n_g0, O(0.1 * scale) for scale. iminuit picks reasonable
        # defaults but explicit hints help on this multi-scale problem.
        m.errors['e_j_hz'] = max(abs(ej0) * 0.05, 1e7)
        m.errors['e_c_hz'] = max(abs(ec0) * 0.05, 1e6)
        m.errors['n_g0'] = 0.01
        m.errors['scale'] = max(abs(scale0) * 0.1, 1e-3)
        m.errors['baseline'] = max(abs(baseline0) * 0.1,
                                   np.std(cq_measured) * 1e-2, 1e-9)

        m.limits['e_j_hz'] = e_j_bounds_hz
        m.limits['e_c_hz'] = e_c_bounds_hz
        m.limits['n_g0'] = n_g0_bounds
        m.limits['scale'] = scale_bounds
        m.limits['baseline'] = baseline_bounds

        m.fixed['n_g0'] = not fit_offset
        m.fixed['scale'] = not fit_scale
        m.fixed['baseline'] = not fit_baseline
        if not fit_offset:
            m.values['n_g0'] = 0.0
        if not fit_scale:
            m.values['scale'] = fixed_scale
        if not fit_baseline:
            m.values['baseline'] = fixed_baseline

        # Coarse 2D pre-scan in (ratio, E_C) to locate the chi² basin
        # before iminuit refines locally. The C_Q(n_g) chi² landscape
        # has a shallow valley along the (E_J, E_C) "weak shape"
        # direction; without this scan iminuit can settle in a
        # different basin from the global minimum depending on noise.
        ratio_grid = np.geomspace(2.0, 50.0, 11)
        ec_grid = np.geomspace(max(e_c_bounds_hz[0], ec0 / 3),
                               min(e_c_bounds_hz[1], ec0 * 3), 7)
        best_fval = np.inf
        best_start = (ej0, ec0, ng0_0, scale0, baseline0)
        for ec_try in ec_grid:
            for ratio in ratio_grid:
                ej_try = ec_try * ratio
                if not (e_j_bounds_hz[0] <= ej_try <= e_j_bounds_hz[1]):
                    continue
                f = cost(ej_try, ec_try, ng0_0, scale0, baseline0)
                if f < best_fval:
                    best_fval = f
                    best_start = (ej_try, ec_try, ng0_0,
                                  scale0, baseline0)
        m.values['e_j_hz'] = float(np.clip(best_start[0], *e_j_bounds_hz))
        m.values['e_c_hz'] = float(np.clip(best_start[1], *e_c_bounds_hz))
        m.values['n_g0'] = ng0_0
        if fit_scale:
            m.values['scale'] = scale0
        if fit_baseline:
            m.values['baseline'] = baseline0

        # MIGRAD + HESSE refinement.
        m.migrad()
        m.hesse()

        e_j = float(m.values['e_j_hz'])
        e_c = float(m.values['e_c_hz'])
        n_g0 = float(m.values['n_g0'])
        scale = float(m.values['scale'])
        baseline = float(m.values['baseline'])

        cq_fit = model_eval(e_j, e_c, n_g0, scale, baseline)

        sig_ej = float(m.errors['e_j_hz'])
        sig_ec = float(m.errors['e_c_hz'])

        # Propagate the (E_J, E_C) HESSE covariance to E_J/E_C:
        # σ_r² = (σ_J/E_C)² + (E_J σ_C / E_C²)² − 2 (E_J/E_C³) cov_JC.
        try:
            cov_jc = float(m.covariance['e_j_hz', 'e_c_hz'])
        except Exception:
            cov_jc = 0.0
        ratio_var = ((sig_ej / e_c) ** 2
                     + (e_j * sig_ec / e_c ** 2) ** 2
                     - 2.0 * (e_j / e_c ** 3) * cov_jc)
        sig_ratio = float(np.sqrt(max(ratio_var, 0.0)))

        errs = {
            'e_j_hz': sig_ej,
            'e_c_hz': sig_ec,
            'n_g0': float(m.errors['n_g0']) if fit_offset else 0.0,
            'scale': float(m.errors['scale']) if fit_scale else 0.0,
            'baseline': (float(m.errors['baseline'])
                         if fit_baseline else 0.0),
            'ej_ec_ratio': sig_ratio,
        }

        return {
            'e_j_hz': e_j,
            'e_c_hz': e_c,
            'n_g0': n_g0,
            'scale': scale,
            'baseline': baseline,
            'ej_ec_ratio': e_j / e_c,
            'errors': errs,
            'cq_fit': cq_fit,
            'residuals': cq_measured - cq_fit,
            'fval': float(m.fval),
            'minuit': m,
        }

    def plot_quantum_capacitance(self, offset_charges=None, c_g_f=None,
                                 charge_cutoff=18, figsize=(4, 3)):
        """
        Plot quantum capacitance C_Q(n_g) for both parities.

        Parameters
        ----------
        offset_charges : array_like, optional
            Offset-charge grid; defaults to linspace(0, 1, 500).
        c_g_f : float, optional
            Gate capacitance [F]. If given, y-axis is in fF; if None
            the intrinsic ∂²E/∂n_g² is plotted in GHz.
        charge_cutoff : int, optional
            Charge-basis cutoff, default 18.
        figsize : tuple, optional
            Figure size, default (4, 3).

        Returns
        -------
        fig, ax : matplotlib figure and axes
        """
        if offset_charges is None:
            offset_charges = np.linspace(0, 1, 500)

        cq_even, cq_odd = self.compute_quantum_capacitance(
            offset_charges, c_g_f=c_g_f, charge_cutoff=charge_cutoff,
        )

        if c_g_f is None:
            y_even = cq_even / 1e9
            y_odd = cq_odd / 1e9
            ylabel = r'$-\partial^2 E / \partial n_g^2$ [GHz]'
        else:
            y_even = cq_even * 1e15
            y_odd = cq_odd * 1e15
            ylabel = r'$C_Q$ [fF]'

        with plt.style.context(self._style_path):
            fig, ax = plt.subplots(figsize=figsize)
            ax.plot(offset_charges, y_even, color='tab:red',
                    linewidth=2, label='even')
            ax.plot(offset_charges, y_odd, color='tab:blue',
                    linewidth=2, label='odd')
            ax.set_xlabel(r'Offset Charge [$C_g V_g / 2e$]')
            ax.set_ylabel(ylabel)
            ax.set_xlim([offset_charges[0], offset_charges[-1]])
            title_parts = [
                f'$\\xi={self.ej_ec_ratio:.1f}$',
                f'$E_J={self.e_j_hz/1e9:.2f}$ GHz',
                f'$E_C={self.e_c_hz/1e9:.3f}$ GHz',
            ]
            ax.set_title(', '.join(title_parts), fontsize=7)
            ax.minorticks_on()
            ax.legend(loc='best', fontsize=7)
            ax.grid(alpha=0.3)

        return fig, ax

    def plot_capacitance_fit(self, offset_charges, cq_measured,
                             fit_result, units='Hz', figsize=(4, 4)):
        """
        Plot a measured C_Q trace with the best-fit overlay and residuals.

        Parameters
        ----------
        offset_charges : array_like
            n_g values at which the data was sampled.
        cq_measured : array_like
            Measured C_Q values (same units as `fit_result['cq_fit']`).
        fit_result : dict
            Output of `fit_quantum_capacitance`.
        units : str, optional
            Unit string used to label the y-axis of both panels. Default
            'Hz' (intrinsic curvature). Use e.g. 'F', 'fF', or 'aF' when
            the input was in absolute capacitance units.
        figsize : tuple, optional
            Figure size, default (4, 4).

        Returns
        -------
        fig, (ax_top, ax_bot) : figure and (data, residual) axes.
        """
        offset_charges = np.asarray(offset_charges, dtype=float)
        cq_measured = np.asarray(cq_measured, dtype=float)
        cq_fit = fit_result['cq_fit']
        resid = fit_result['residuals']

        with plt.style.context(self._style_path):
            fig, (ax_top, ax_bot) = plt.subplots(
                2, 1, figsize=figsize, sharex=True,
                gridspec_kw={'height_ratios': [3, 1]},
            )
            ax_top.plot(offset_charges, cq_measured, 'o',
                        markersize=3, color='black', label='measured')
            ax_top.plot(offset_charges, cq_fit, '-',
                        color='tab:red', linewidth=2, label='fit')
            ax_top.set_ylabel(rf'$C_Q$ [{units}]')
            errs = fit_result.get('errors', {}) or {}
            sig_ratio = errs.get('ej_ec_ratio', 0.0)
            sig_ej = errs.get('e_j_hz', 0.0)
            sig_ec = errs.get('e_c_hz', 0.0)
            sig_ng0 = errs.get('n_g0', 0.0)
            sig_bl = errs.get('baseline', 0.0)
            baseline = fit_result.get('baseline', 0.0)
            title_parts = [
                rf"$E_J/E_C={fit_result['ej_ec_ratio']:.2f}"
                rf"\pm{sig_ratio:.2f}$",
                rf"$E_J=({fit_result['e_j_hz']/1e9:.3f}"
                rf"\pm{sig_ej/1e9:.3f})$ GHz",
                rf"$E_C=({fit_result['e_c_hz']/1e9:.4f}"
                rf"\pm{sig_ec/1e9:.4f})$ GHz",
                rf"$n_{{g0}}={fit_result['n_g0']:+.3f}"
                rf"\pm{sig_ng0:.3f}$",
            ]
            # Only annotate baseline when it was actually fit (nonzero
            # uncertainty) or when the user supplied a nonzero fixed
            # value — otherwise it just clutters the title.
            if sig_bl > 0 or baseline != 0:
                title_parts.append(
                    rf"$y_0=({baseline:+.3g}\pm{sig_bl:.2g})$ {units}"
                )
            ax_top.set_title(', '.join(title_parts), fontsize=7)
            ax_top.legend(loc='best', fontsize=7)
            ax_top.minorticks_on()
            ax_top.grid(alpha=0.3)

            ax_bot.plot(offset_charges, resid, '.',
                        markersize=3, color='tab:gray')
            ax_bot.axhline(0, color='black', linewidth=0.8, alpha=0.6)
            ax_bot.set_xlabel(r'Offset Charge [$C_g V_g / 2e$]')
            ax_bot.set_ylabel(f'residual [{units}]')
            ax_bot.minorticks_on()
            ax_bot.grid(alpha=0.3)

        return fig, (ax_top, ax_bot)

    def plot_likelihood_landscape(self, offset_charges, cq_measured,
                                  fit_result, parity='odd',
                                  n_grid=21, span_factor=1.5,
                                  charge_cutoff=10, cq_sigma=None,
                                  figsize=(4.5, 4)):
        """
        Plot the chi² landscape in the (E_J, E_C) plane around a fit.

        Visualises the likelihood-function curvature so the user can
        see the shallow valley that limits joint (E_J, E_C) recovery
        from a C_Q(n_g) fit. n_g0 and `scale` are held fixed at their
        fitted values; only E_J and E_C are scanned.

        Parameters
        ----------
        offset_charges : array_like
            Same n_g grid that was used in `fit_quantum_capacitance`.
        cq_measured : array_like
            The measured trace.
        fit_result : dict
            Output of `fit_quantum_capacitance`. Used as the centre of
            the scan and for n_g0 / scale.
        parity : {'odd', 'even'}, optional
            Parity branch the data was fit on. Must match the fit.
        n_grid : int, optional
            Grid resolution per axis. Default 21 (→ 441 evaluations).
        span_factor : float, optional
            (E_J, E_C) range as multiplicative factor around the fit.
            E.g. 1.5 means [fit/1.5, fit*1.5] log-spaced. Default 1.5.
        charge_cutoff : int, optional
            Charge-basis cutoff for the inner solves. Default 10
            (faster than the fit's 18; sufficient for visualisation
            in the QPD/transmon regime).
        cq_sigma : float or array_like, optional
            Per-point measurement uncertainty used to normalise the
            chi². Default uses the same heuristic as
            `fit_quantum_capacitance`.
        figsize : tuple, optional
            Figure size, default (4.5, 4).

        Returns
        -------
        fig, ax : matplotlib figure and axes.

        Notes
        -----
        Confidence-region contours overlay the chi² surface at
        Δχ² = {2.30, 6.18, 11.83}, which correspond to 68 %, 95 %,
        and 99 % joint coverage in 2 parameters.
        """
        offset_charges = np.asarray(offset_charges, dtype=float)
        cq_measured = np.asarray(cq_measured, dtype=float)

        if parity not in ('odd', 'even'):
            raise ValueError(f"Invalid parity: {parity!r}")

        if cq_sigma is None:
            sigma = max(np.std(cq_measured) * 1e-2,
                        np.max(np.abs(cq_measured)) * 1e-3)
        else:
            sigma = float(np.asarray(cq_sigma).mean())

        ej_fit = float(fit_result['e_j_hz'])
        ec_fit = float(fit_result['e_c_hz'])
        ng0_fit = float(fit_result['n_g0'])
        scale_fit = float(fit_result['scale'])
        baseline_fit = float(fit_result.get('baseline', 0.0))

        ej_grid = np.geomspace(ej_fit / span_factor,
                               ej_fit * span_factor, n_grid)
        ec_grid = np.geomspace(ec_fit / span_factor,
                               ec_fit * span_factor, n_grid)

        chi2 = np.empty((n_grid, n_grid))
        for i, ec in enumerate(ec_grid):
            for j, ej in enumerate(ej_grid):
                tmp = QPD.__new__(QPD)
                tmp.e_j_hz = ej
                tmp.e_c_hz = ec
                tmp.e_j_ev = ej * QPD.PLANCK_EV_S
                tmp.e_c_ev = ec * QPD.PLANCK_EV_S
                tmp.delta_l_ev = 0.0
                tmp.delta_r_ev = 0.0
                cq_e, cq_o = QPD.compute_quantum_capacitance(
                    tmp, offset_charges - ng0_fit, c_g_f=None,
                    charge_cutoff=charge_cutoff,
                )
                model = (scale_fit * (cq_o if parity == 'odd' else cq_e)
                         + baseline_fit)
                resid = (model - cq_measured) / sigma
                chi2[i, j] = float(np.sum(resid * resid))

        chi2_min = chi2.min()
        delta_chi2 = chi2 - chi2_min

        # 2D joint-coverage levels
        levels = (2.30, 6.18, 11.83)
        labels = {2.30: r'$1\sigma$', 6.18: r'$2\sigma$',
                  11.83: r'$3\sigma$'}

        # Visualise log10(1 + Δχ²) so both the deep valley *and* the
        # walls remain readable instead of saturating to one colour.
        log_dchi2 = np.log10(1.0 + delta_chi2)

        with plt.style.context(self._style_path):
            fig, ax = plt.subplots(figsize=figsize)
            EJ, EC = np.meshgrid(ej_grid / 1e9, ec_grid / 1e9)

            pcm = ax.pcolormesh(EJ, EC, log_dchi2,
                                cmap='magma_r', shading='auto')
            cbar = fig.colorbar(pcm, ax=ax, pad=0.02)
            cbar.set_label(r'$\log_{10}(1+\Delta\chi^2)$', fontsize=8)

            cs = ax.contour(EJ, EC, delta_chi2,
                            levels=levels, colors='white',
                            linewidths=1.0)
            ax.clabel(cs, inline=True, fontsize=6,
                      fmt={lvl: labels[lvl] for lvl in levels})

            # Constant-ratio reference line through the fit point
            ratio = ej_fit / ec_fit
            ec_line = np.array([ec_grid[0], ec_grid[-1]])
            ax.plot(ratio * ec_line / 1e9, ec_line / 1e9,
                    '--', color='tab:cyan', linewidth=0.8,
                    alpha=0.7,
                    label=f'$E_J/E_C={ratio:.2f}$')

            # Mark the fit point on top
            ax.plot(ej_fit / 1e9, ec_fit / 1e9, '+',
                    color='tab:cyan', markersize=12,
                    markeredgewidth=2.0, label='fit')

            ax.set_xlabel(r'$E_J$ [GHz]')
            ax.set_ylabel(r'$E_C$ [GHz]')
            ax.set_xscale('log'); ax.set_yscale('log')
            ax.set_title(
                rf'$\chi^2$ landscape; $\chi^2_\text{{min}}$ = '
                rf'{chi2_min:.1f}', fontsize=8)
            leg = ax.legend(loc='lower right', fontsize=7,
                            facecolor='white', edgecolor='black',
                            framealpha=0.9, labelcolor='black')
            leg.get_frame().set_linewidth(0.5)
            ax.minorticks_on()

        return fig, ax

    def plot_all(self, offset_charges=None, coupling_g_hz=150e6,
                readout_freq_hz=7.0e9, num_levels=5):
        """
        Generate all standard plots for QPD transmon analysis
        
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
        print(f"QPD Transmon Parameters:")
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
    """Example usage of QPD class"""
    # Example from the MATLAB script (WashU parameters)
    ej_ec_ratio = 12
    e_j_hz = 8.335e9  # ~8.3 GHz (0.4 K·kB)
    e_c_hz = e_j_hz / ej_ec_ratio  # Hz
    
    # Create QPD instance
    qpd = QPD(
        e_j_hz=e_j_hz,
        e_c_hz=e_c_hz,
        temperature_k=0.012,  # 12 mK
        r_n_ohm=27e3  # 27 kΩ
    )
    
    # Generate all plots
    figs = qpd.plot_all(
        coupling_g_hz=150e6,  # 150 MHz
        readout_freq_hz=7.0e9,  # 7 GHz
        num_levels=5
    )
    
    plt.show()
    
    return qpd, figs


if __name__ == "__main__":
    qpd, figs = main()

