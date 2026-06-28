# ASSERT_CONVENTION: units_policy=Hz_api_radps_solver, two_pi_site=single_at_hz_to_qobj_wrap, input_output=gardiner_collett_minus_i_kappa_c_port_total_kappa_dissipator, drive=eps_a_plus_adag_nonRWA, gksl=L_rho_Ldag_minus_half_anticommutator, tensor_order=qubit_major_kron_Hq_Ir, coupling_form=g_nhat_minus_ng_a_plus_adag_nonRWA_Delta_omegaq_minus_omegar, dissipator_basis=dressed_basis_C5_undriven_eigenbasis_cavity_bare_total_kappa, frame_policy=lab_frame_static_nonRWA_coupling_floquet_steadystate_fourier_rho1_sideband
"""Driven-dissipative open-system circuit-QED solver foundation (Phase 1, Plan 01).

This module holds ONLY the foundation pieces proven in the SIMU-01 cavity unit gate:

  * the SINGLE 2*pi conversion site at the Hz -> Qobj wrap (C1),
  * a cavity-only (no qubit) driven-damped steady-state Liouvillian builder,
  * the Gardiner-Collett input-output S21 map (C2).

The qubit Hamiltonian wrap, dressed-basis dissipators (C5), and the PR #19 (VALD-01)
comparison are deliberately NOT here -- they are Plan 02 and Plan 03.

Convention lock (CONVENTIONS C1-C3, the registry SIMU-01 enforces):

  C1 Units.   All device quantities (E_J, E_C, f_r, g, kappa, kappa_c, Gamma, eps)
              are ordinary frequency in Hz (= H/h) at the API/storage boundary.
              Inside QuTiP everything is angular frequency (rad/s, hbar=1). The
              2*pi conversion happens at EXACTLY ONE site -- ``hz_to_qobj`` -- and
              the collapse-operator rate carries the 2*pi (``sqrt(2*pi*rate_hz)``).

  C2 S21.     S21(omega) = 1 - i*sqrt(kappa_c)*<a>/alpha_in (Gardiner-Collett, -i
              drive phase), a_out = a_in - i*sqrt(kappa_c)*<a>. The damping
              dissipator uses the TOTAL kappa; the input-output port/prefactor uses
              kappa_c. Phase-1 default is kappa_c = kappa, but the two are kept as
              SEPARATE named variables so the distinction never silently collapses.

  C3 Drive.   H_d = eps*(a + a_dag), non-RWA cavity-port drive. eps = (kappa/2)*sqrt(<n>),
              alpha_in = eps/sqrt(kappa_c). On resonance with eps = kappa/2 the empty
              cavity has <n> = 1.

Closed forms reproduced by the cavity-only solve (the SIMU-01 asserts in
``the cavity unit-gate checks``):

    <a>(omega_d) = -i*eps / (i*(omega_r - omega_d) + kappa/2)        # on res: -2i*eps/kappa
    <n>(omega_d) =  eps**2 / ((omega_r - omega_d)**2 + (kappa/2)**2) # FWHM (in omega_d) = kappa
    S21(omega_r) = 1 - 2*kappa_c/kappa                               # = -1 when kappa_c = kappa

Reference: SQUAT App. E (Magoon et al., arXiv:2601.16261), adapted so the port field
is the cavity ``a`` and the dissipator uses total kappa (C2); CONVENTIONS C1-C3.
"""

from __future__ import annotations

import numpy as np
import scipy.linalg as _sla
from qutip import (
    Qobj,
    destroy,
    expect,
    liouvillian,
    operator_to_vector,
    qeye,
    qeye_like,
    qzero_like,
    steadystate,
    tensor,
    vector_to_operator,
)
from qutip.core import data as _data

__all__ = [
    "hz_to_qobj",
    "collapse",
    "cavity_only_liouvillian",
    "s21",
    "wrap_cqed_hamiltonian",
    "cavity_operator",
    "dressed_jump_operators",
    "dressed_jminus_roundtrip_check",
    "displacement_alphas",
    "charge_drive_operator",
    "driven_steadystate",
    "nodrive_precondition",
    "driven_s21_sweep",
    "extract_peak_hz",
    "FloquetSidebandSolver",
    "vald01_peak_scan",
    "vald01_gate_point",
    "vald01_measured_tolerance",
    "vald01_run_gate",
    "vald01_eps_extrapolation",
    "mesolve_branch_run",
    "semiclassical_floquet_onset",
    "sambe_floquet_quasienergies",
    "sambe_floquet_branch_scan",
    "cavity_qfunc",
    "pointer_snr",
]


def hz_to_qobj(H_hz, dims=None) -> Qobj:
    """The SINGLE 2*pi conversion site (C1): Hz -> rad/s Qobj.

    Every conversion of a frequency-dimensioned Hamiltonian from the Hz API regime
    into the rad/s solver regime flows through here and nowhere else. ``2*np.pi``
    multiplies a frequency-carrying operator at this one site.

    Parameters
    ----------
    H_hz : numpy.ndarray or Qobj
        Hamiltonian (or any frequency-dimensioned operator) in Hz (= H/h).
    dims : list, optional
        QuTiP composite-space dims, e.g. ``[[n_qubit, n_fock], [n_qubit, n_fock]]``.
        Ignored if ``H_hz`` is already a Qobj (its dims are kept).

    Returns
    -------
    Qobj or float
        ``2*pi * Qobj(H_hz)`` in rad/s.  If ``H_hz`` is a plain scalar
        (a frequency-dimensioned coefficient such as a detuning or drive
        amplitude that multiplies an operator built elsewhere), returns the
        scalar ``2*pi * H_hz`` -- still the single conversion site, applied to a
        coefficient rather than a full operator matrix.

    Notes
    -----
    Accepting a scalar lets callers write ``hz_to_qobj(eps_hz) * (a + a.dag())``
    and keep the 2*pi factor at this one site (C1) even when the operator
    (here ``a + a.dag()``) is assembled in the joint QuTiP space rather than
    from a raw Hz matrix.
    """
    if isinstance(H_hz, Qobj):
        return 2.0 * np.pi * H_hz
    if np.isscalar(H_hz):
        return 2.0 * np.pi * H_hz
    return 2.0 * np.pi * Qobj(H_hz, dims=dims)


def collapse(rate_hz: float, op: Qobj) -> Qobj:
    """GKSL collapse operator with the rate carrying the 2*pi (C1).

    Returns ``sqrt(2*pi*rate_hz) * op``. With the GKSL dissipator
    ``D[L]rho = L rho L_dag - 0.5{L_dag L, rho}`` this puts a rate ``2*pi*rate_hz``
    (rad/s) into the master equation -- consistent with the single-2*pi-site rule
    (the rate, not the operator, carries the 2*pi).

    Parameters
    ----------
    rate_hz : float
        Decay/dephasing rate in Hz (e.g. the cavity total kappa).
    op : Qobj
        The bare jump operator (e.g. the cavity annihilation ``a``).
    """
    return np.sqrt(2.0 * np.pi * rate_hz) * op


def cavity_only_liouvillian(
    omega_r_hz: float,
    omega_d_hz: float,
    eps_hz: float,
    kappa_hz: float,
    n_fock: int = 30,
):
    """Build and solve the cavity-only driven-damped steady state.

    Rotating frame at the drive frequency omega_d (C-locked): the bare resonator term
    ``omega_r * a_dag a`` becomes ``(omega_r - omega_d) * a_dag a`` and the drive
    ``eps*(a + a_dag)`` is static (the drive's fast 2*omega_d counter-rotating piece is
    dropped at low power; see RESEARCH.md). No qubit, no dressed dissipators -- this is
    the SIMU-01 Checks 1 & 2 system that isolates the 2*pi and sign/sqrt(kappa) bugs.

    .. warning:: (Phase 1.1) This rotating-frame treatment is exact ONLY for this
       cavity-only system (every operator here commutes into the frame rotation;
       the dropped 2*omega_d drive piece is negligible at weak drive). Applying
       the same "rotating frame + static operators" shortcut to the JOINT
       qubit-cavity system while keeping the non-RWA coupling
       ``g(n_hat-n_g)(a+a_dag)`` static is the INDICTED Phase-1 VALD-01 bug
       (fp-rf-shortcut): it promotes ``g*n_hat*a_dag`` to a spurious static
       cavity source (no-drive <n>~1.53, sign-flipped peak). The joint driven
       solver below therefore works in the LAB frame (Floquet periodic steady
       state); this cavity-only path is kept solely for the SIMU-01 regression.

    All inputs are in Hz. The single 2*pi site is ``hz_to_qobj`` (applied to both the
    detuning term and the drive); the collapse rate carries 2*pi via ``collapse``.

    Parameters
    ----------
    omega_r_hz : float
        Bare cavity frequency f_r in Hz.
    omega_d_hz : float
        Drive/probe frequency omega_d in Hz.
    eps_hz : float
        Drive amplitude epsilon in Hz (C3: H_d = eps*(a + a_dag)).
    kappa_hz : float
        TOTAL cavity linewidth kappa in Hz (the damping dissipator rate, C2).
    n_fock : int, optional
        Fock-space cutoff. Default 30 (generous; <n> ~ O(1) at the SIMU-01 drive
        eps = kappa/2 needs only ~10, but a wide cutoff de-risks truncation leakage).

    Returns
    -------
    (rho_ss, a)
        ``rho_ss`` the steady-state density matrix (Qobj), and ``a`` the cavity
        annihilation operator (Qobj) so the caller can take expectations such as
        ``expect(a, rho_ss)`` and ``expect(a.dag() * a, rho_ss)``.
    """
    a = destroy(n_fock)

    # Rotating-frame detuning Hamiltonian in Hz, then wrapped at the single 2*pi site.
    H_detune_hz = (omega_r_hz - omega_d_hz) * (a.dag() * a)
    H = hz_to_qobj(H_detune_hz)

    # Static non-RWA drive (C3), 2*pi applied at the same wrap site.
    H = H + hz_to_qobj(eps_hz * (a + a.dag()))

    # Cavity damping uses the TOTAL kappa (C2); rate carries the 2*pi (C1).
    c_ops = [collapse(kappa_hz, a)]

    rho_ss = steadystate(H, c_ops, method="direct")
    return rho_ss, a


def s21(a_expect: complex, eps_hz: float, kappa_c_hz: float) -> complex:
    """Gardiner-Collett input-output transmission map (C2).

    S21 = 1 - i*sqrt(kappa_c)*<a>/alpha_in, with alpha_in = eps/sqrt(kappa_c).
    The port/prefactor uses ``kappa_c_hz`` (NOT the total kappa used by the
    dissipator); the two are kept as separate named quantities (Phase-1 sets them
    equal). The 2*pi-invariance of S21 (a ratio) means this map needs no 2*pi factor.

    On resonance with <a> = -2i*eps/kappa and kappa_c = kappa this returns
    S21(omega_r) = 1 - 2*kappa_c/kappa = -1 (perfect notch) -- the end-to-end SIMU-01
    gate that catches any sign / 2*pi / sqrt(kappa_c) / -i-phase error at once.

    Parameters
    ----------
    a_expect : complex
        Steady-state intracavity field <a> = Tr(a rho_ss).
    eps_hz : float
        Drive amplitude epsilon in Hz (sets alpha_in = eps/sqrt(kappa_c)).
    kappa_c_hz : float
        Cavity <-> feedline coupling rate kappa_c in Hz (the PORT rate, C2).
    """
    alpha_in = eps_hz / np.sqrt(kappa_c_hz)
    return 1.0 - 1j * np.sqrt(kappa_c_hz) * a_expect / alpha_in


# ===========================================================================
# Plan 02: joint qubit-resonator wrap (C1, C4) and dressed-basis dissipators (C5)
# ===========================================================================
#
# C4  Coupling/detuning.  The joint Hamiltonian is the repo's
#     ``build_cqed_hamiltonian`` (full-cosine multilevel, non-RWA, regime-
#     agnostic):  H/h = diag(omega_q_i) (x) I_r + omega_r I_q (x) a_dag a
#                       + g (n_hat - n_g)_q (x) (a + a_dag).
#     Detuning Delta = omega_q - omega_r (Blais sign); Delta < 0 at the WashU
#     point.  We REUSE that builder verbatim (no charge-basis rebuild,
#     no two-level/fixed-E_J/E_C reduction -- forbidden proxy fp-twolevel-wrap),
#     wrapping its Hz matrix as a Qobj at the single 2*pi site.
#
#     IMPORTANT (C4 matrix element).  The coupling multiplies the FULL operator
#     (n_hat - n_g), whose qubit 0<->1 matrix element n_01 = <0|n_hat|1> is NOT
#     1.  At the WashU QPD point n_01 ~ 0.697.  The physical vacuum-Rabi (Jaynes-
#     Cummings) coupling of the 0<->1 transition is therefore
#         g_01 = g * |n_01|,
#     and the dispersive/Purcell scale is set by g_01, not by the bare device g.
#     The Purcell check (checks/simu01_purcell.py) compares against
#     (g_01/Delta)^2 * kappa = (g |n_01| / Delta)^2 * kappa.
#
# C1  Tensor order.  ``build_cqed_hamiltonian`` uses np.kron(H_qubit, I_r), i.e.
#     QUBIT-MAJOR with joint index = i*(n_photon+1) + n.  In QuTiP this is
#     tensor(A_qubit, B_resonator), dims = [[n_qubit, n_fock], [n_qubit, n_fock]].
#     The bare cavity field is a = tensor(qeye(n_qubit), destroy(n_fock)).
#
# C5  Dissipators DRESSED-FROM-THE-START.  Qubit jump operators are built in the
#     FULL coupled UNDRIVEN joint eigenbasis (solve_cqed_eigensystem), not the
#     bare basis (forbidden proxy fp-bare-dissipator -- a bare dissipator gives
#     Gamma_Purcell = 0).  The cavity damping stays kappa*D[a] in the BARE field
#     a with the TOTAL kappa.  GKSL form D[L]rho = L rho L_dag - 0.5{L_dag L, rho};
#     the rate carries the 2*pi via ``collapse`` (C1).
# ===========================================================================


def wrap_cqed_hamiltonian(
    qpd,
    offset_charge: float,
    coupling_g_hz: float,
    readout_freq_hz: float,
    parity: str = "odd",
    n_qubit: int = 16,
    n_fock: int = 13,
    charge_cutoff: int = 30,
):
    """Wrap ``build_cqed_hamiltonian`` (Hz) as a rad/s QuTiP ``Qobj`` (C1, C4).

    REUSES the repo's validated joint Hamiltonian builder verbatim -- no
    charge-basis rebuild, no two-level / fixed-E_J/E_C reduction (forbidden
    proxy ``fp-twolevel-wrap``).  The returned NumPy Hz matrix is wrapped at the
    single 2*pi site (``hz_to_qobj``) with the qubit-major composite dims that
    match ``np.kron(H_qubit, I_r)`` (joint index = i*n_fock + n).

    Parameters
    ----------
    qpd : qpd.theory.transmon.QPD
        Transmon instance carrying E_J, E_C (full-cosine CPB).
    offset_charge : float
        Dimensionless offset charge n_g (period 1 = 2e; parity adds 0.5 for
        'odd', C6 -- already wired inside build_cqed_hamiltonian).
    coupling_g_hz : float
        Bare device coupling g [Hz].  NOTE: the qubit 0<->1 vacuum-Rabi coupling
        is g_01 = g*|n_01| with n_01 the qubit number matrix element (~0.697 at
        WashU), not g itself.
    readout_freq_hz : float
        Bare resonator frequency omega_r [Hz].
    parity : {'odd', 'even'}
        CPB parity (odd shifts n_g by +0.5, C6).
    n_qubit : int
        Number of transmon eigenstates kept.
    n_fock : int
        Fock dimension (= n_photon + 1).
    charge_cutoff : int
        CPB charge-basis cutoff passed through to the builder.

    Returns
    -------
    H : Qobj
        Joint Hamiltonian in rad/s, dims = [[n_qubit, n_fock], [n_qubit, n_fock]].
    components : dict
        ``{'H_hz', 'qubit_freqs_hz', 'n_shift_mat', 'basis_labels',
        'n_qubit', 'n_fock'}`` from ``return_components`` -- carried so the
        caller can build the dressed jump operators (Task 2) without re-solving.
    """
    n_photon = n_fock - 1
    H_hz, qubit_freqs_hz, n_shift_mat, basis_labels = qpd.build_cqed_hamiltonian(
        offset_charge,
        coupling_g_hz,
        readout_freq_hz,
        parity=parity,
        n_qubit=n_qubit,
        n_photon=n_photon,
        rwa=False,
        charge_cutoff=charge_cutoff,
        return_components=True,
    )
    dims = [[n_qubit, n_fock], [n_qubit, n_fock]]
    # Single 2*pi site (C1): Hz -> rad/s, qubit-major composite dims.
    H = hz_to_qobj(np.asarray(H_hz, dtype=float), dims=dims)
    components = {
        "H_hz": np.asarray(H_hz, dtype=float),
        "qubit_freqs_hz": np.asarray(qubit_freqs_hz, dtype=float),
        "n_shift_mat": np.asarray(n_shift_mat, dtype=float),
        "basis_labels": basis_labels,
        "n_qubit": n_qubit,
        "n_fock": n_fock,
    }
    return H, components


def cavity_operator(n_qubit: int, n_fock: int) -> Qobj:
    """Bare cavity annihilation ``a`` in the joint space, qubit-major (C1).

    ``a = tensor(qeye(n_qubit), destroy(n_fock))`` -- consistent with
    ``build_cqed_hamiltonian``'s ``np.kron(H_qubit, I_r)`` ordering so that the
    joint index is i*n_fock + n.  This is the BARE field that the cavity
    dissipator kappa*D[a] acts on (C5: cavity damping stays in the bare field).
    """
    return tensor(qeye(n_qubit), destroy(n_fock))


def dressed_jump_operators(
    qpd,
    offset_charge: float,
    coupling_g_hz: float,
    readout_freq_hz: float,
    parity: str = "odd",
    n_qubit: int = 16,
    n_fock: int = 13,
    charge_cutoff: int = 30,
):
    """C5 dressed-basis qubit jump operators from the undriven joint eigenbasis.

    Diagonalize the FULL coupled UNDRIVEN qubit+resonator Hamiltonian via the
    repo's ``solve_cqed_eigensystem``, then build the qubit relaxation and
    dephasing operators by the *projected-bare-operator* form (RESEARCH.md
    recommendation, captures Purcell automatically):

        J_-^dressed   = sum_{E_k < E_l} <psi_k|(sigma_-^q (x) I_r)|psi_l> |psi_k><psi_l|
        Sz_dressed    = U^dag (sigma_z^q (x) I_r) U   (kept full; ~diagonal in
                        the dressed basis to leading order)

    where ``U`` has the dressed eigenvectors as columns and the energy-LOWERING
    restriction (E_k < E_l) on J_- keeps only relaxation terms.  Building these
    in the *dressed* basis (not the bare basis -- forbidden proxy
    ``fp-bare-dissipator``) is what lets cavity damping kappa*D[a] leak qubit
    decay through the cavity (the Purcell effect): the dressed qubit-excited
    state carries a small bare-photon admixture ~ g*n_01/Delta.

    Bare qubit operators
    --------------------
    The qubit block is the *energy eigenbasis* (states ordered by energy), so
    the bare qubit lowering operator is the unit lower-shift
    ``sigma_-^q = sum_i |i-1><i|`` (each |i> -> |i-1>), and
    ``sigma_z^q = diag(+1, -1, -1, ..., -1)`` projecting the 0<->1 qubit
    pseudospin (ground +1, all excited -1; C5 Gamma_phi/2 D[sigma_z]).  These
    are the standard pseudospin operators; the physical 0<->1 coupling strength
    (the g*n_01 matrix element) lives in the HAMILTONIAN, so the Purcell rate is
    set by the dressed admixture regardless of the jump-operator normalization.

    Returns
    -------
    ops : dict
        ``{'J_minus_dressed', 'Sz_dressed', 'J_minus_bare', 'Sz_bare', 'a',
        'evals_hz', 'evecs', 'labels', 'overlaps', 'n_qubit', 'n_fock'}``.
        ``J_minus_dressed`` / ``Sz_dressed`` are the C5 collapse-operator
        FACTORS (multiply by ``collapse(rate_hz, op)`` to form the GKSL jump
        operator); ``*_bare`` are the bare-basis counterparts kept for the
        negative-control comparison (they give Gamma_Purcell = 0).
    """
    n_photon = n_fock - 1
    evals_hz, evecs, labels, overlaps = qpd.solve_cqed_eigensystem(
        offset_charge,
        coupling_g_hz,
        readout_freq_hz,
        parity=parity,
        n_qubit=n_qubit,
        n_photon=n_photon,
        rwa=False,
        charge_cutoff=charge_cutoff,
    )
    dims = [[n_qubit, n_fock], [n_qubit, n_fock]]

    # --- Bare qubit pseudospin operators (energy-ordered eigenbasis) ---------
    # sigma_-^q : the energy-LOWERING qubit shift |i> -> |i-1|.  As a matrix
    # acting on column vectors, lowering means (sigma_-)[i-1, i] = 1, i.e. the
    # SUPER-diagonal (k=+1).  [Plan-03 DEVIATION Rule 1 fix: the original
    # k=-1 (sub-diagonal) was a RAISING operator |i> -> |i+1>; combined with the
    # subsequent energy-lowering triu mask it left an almost-null J_- that could
    # NOT relax the qubit, so the driven steady state pinned the qubit in an
    # excited, photon-populated state (P_exc ~ 0.9) instead of the |0,n> ground
    # branch PR #19 reads chi from.  Plan-02's Purcell check never exercised J_-
    # as a relaxation channel (it used kappa*D[a] only), so the sign survived.]
    sigma_minus_q = np.diag(np.ones(n_qubit - 1), k=+1)
    # sigma_z^q : ground +1, all excited -1 (0<->1 pseudospin; C5 dephasing).
    sigma_z_diag = -np.ones(n_qubit)
    sigma_z_diag[0] = 1.0
    sigma_z_q = np.diag(sigma_z_diag)

    I_r = np.eye(n_fock)
    J_minus_bare_full = np.kron(sigma_minus_q, I_r)   # qubit-major (C1)
    Sz_bare_full = np.kron(sigma_z_q, I_r)

    J_minus_bare = Qobj(J_minus_bare_full, dims=dims)
    Sz_bare = Qobj(Sz_bare_full, dims=dims)

    # --- Project into the dressed eigenbasis (C5) ----------------------------
    # U columns = dressed eigenvectors in the bare qubit-major basis.
    U = evecs                                   # (D, D)
    # Bare operator in the dressed basis: O_dressed = U^dag O_bare U.
    Jd_full = U.conj().T @ J_minus_bare_full @ U
    Sz_dressed_full = U.conj().T @ Sz_bare_full @ U

    # Energy-LOWERING restriction on J_-: keep <k|...|l> only for E_k < E_l.
    # evals_hz is sorted ascending (eigh), so row index < col index => E_k < E_l.
    lower_mask = np.zeros((evals_hz.size, evals_hz.size), dtype=bool)
    iu = np.triu_indices(evals_hz.size, k=1)    # k (row) < l (col)
    lower_mask[iu] = True
    Jd_lowering = Jd_full * lower_mask

    # Re-express the dressed-basis operators back in the bare basis so they can
    # be applied directly in the joint QuTiP space (same basis as H and a).
    J_minus_dressed = Qobj(U @ Jd_lowering @ U.conj().T, dims=dims)
    Sz_dressed = Qobj(U @ Sz_dressed_full @ U.conj().T, dims=dims)

    a = cavity_operator(n_qubit, n_fock)

    return {
        "J_minus_dressed": J_minus_dressed,
        "Sz_dressed": Sz_dressed,
        "J_minus_bare": J_minus_bare,
        "Sz_bare": Sz_bare,
        "a": a,
        "evals_hz": evals_hz,
        "evecs": evecs,
        "labels": labels,
        "overlaps": overlaps,
        "n_qubit": n_qubit,
        "n_fock": n_fock,
    }


def dressed_jminus_roundtrip_check(
    qpd,
    offset_charge: float,
    parity: str,
    *,
    coupling_g_hz: float,
    readout_freq_hz: float,
    kappa_hz: float,
    gamma_r_hz: float,
    n_qubit: int = 8,
    n_fock: int = 6,
    charge_cutoff: int = 30,
    pexc_tol: float = 0.01,
    mask_noise_tol: float = 1e-12,
    roundtrip_tol: float = 1e-12,
    lowering_weight_min: float = 0.99,
    assert_pass: bool = True,
):
    """Dressed J_- GKSL round-trip revalidation (Phase 1.1 Plan 03, test-gksl-roundtrip).

    The Plan-03 J_- direction fix (k=+1 super-diagonal, raising -> lowering;
    see the DEVIATION note in :func:`dressed_jump_operators`) is in place and
    the Purcell gate passed, but Purcell exercised ``kappa*D[a]`` ONLY --
    ``D[J_-^dressed]`` as a relaxation channel was never round-trip-checked in
    the bare basis (fp-bare-dissipator: the gap this check closes).  Three
    checks on the DRESSED operator (never a bare-basis substitute):

    1. RELAXATION.  Solve the no-drive steady state under
       ``c_ops = [collapse(kappa, a), collapse(Gamma_r, J_-^dressed)]``
       (the original DERV-01 section-6.2 collapse set -- no dephasing term, so
       J_-^dressed is the ONLY qubit relaxation path) and confirm the
       dressed-qubit excited population

           <P_exc> = sum_{k: qubit label(k) != 0} <psi_k| rho_ss |psi_k>

       is small [dimensionless].  The Phase-1 Plan-03 measurement of the same
       quantity was ~0.007 (in the buggy rotating frame, where the qubit was
       correctly relaxed even while the cavity floor was corrupted at
       <n>~1.53); in the cured lab frame it should be at most that.  With the
       OLD direction bug J_- was nearly null and the qubit pinned excited
       (P_exc ~ 0.9) -- this is the failure mode the check disconfirms.

    2. ENERGY-LOWERING ONLY.  Express the returned ``J_minus_dressed`` back in
       the dressed eigenbasis (Jd = U^dag J U) and verify it has ONLY
       energy-lowering entries: <psi_k|J|psi_l> = 0 for E_k >= E_l (strictly
       lower triangle + diagonal, evals ascending) above numerical noise.  A
       residual energy-raising entry would signal an incomplete direction fix.
       The NON-CIRCULAR companion is the lowering-weight fraction
       ||Jd_masked||_F / ||U^dag J_bare U||_F: with the old k=-1 bug the
       projected bare operator was dominantly RAISING (lower-triangular), so
       the energy-lowering mask annihilated almost all of it (near-null J_-);
       with the fix the mask retains essentially all the weight (the removed
       raising part is only the O(g_01/Delta) dressed admixture).

    3. ROUND-TRIP.  Transform the bare lowering operator into the dressed
       basis and back, ``J_rt = U (U^dag J_bare U) U^dag``, and confirm
       ``J_rt == J_bare`` to numerical tolerance (basis-transform bookkeeping:
       catches any vec-ordering / conjugation / non-unitarity error in the
       dressed projection).  U^dag U = I is checked alongside.

    Inherited validation context (Plans 01-02 / 01.1-01, NOT re-run here):
    Purcell = 35.0 Hz = 1.025*(g_01/Delta)^2*kappa with g_01 = g*|n_01|
    (n_01 = 0.6968) at the QPD4 point -- the rate-level validation of the same
    dressed channel.  This check validates the OPERATOR-level GKSL structure
    that Purcell never exercised.

    Parameters
    ----------
    qpd : qpd.theory.transmon.QPD
        Transmon instance (E_J, E_C in Hz).
    offset_charge, parity
        Gate point (C6: odd shifts n_g by +0.5 inside the builder).
    coupling_g_hz, readout_freq_hz, kappa_hz, gamma_r_hz : float
        Device coupling g, bare resonator f_r, TOTAL cavity kappa, and the
        dressed qubit relaxation rate Gamma_r, all in Hz (C1).
    n_qubit, n_fock, charge_cutoff : int
        Truncation (defaults = the Plan-01/02 gating dims nq=8, nf=6).
    pexc_tol : float
        Bound on <P_exc> (default 0.01, consistent with the ~0.007 Plan-03
        scale; the lab-frame value is expected far below).
    mask_noise_tol : float
        Relative bound on energy-raising / diagonal entries of the dressed
        J_- (numerical-noise scale).
    roundtrip_tol : float
        Relative bound on the bare->dressed->bare round-trip error and on
        ||U^dag U - I||.
    lowering_weight_min : float
        Lower bound on the lowering-weight fraction (near-null J_- detector).
    assert_pass : bool
        If True (default), raise ``ValueError`` on any failed check.

    Returns
    -------
    dict with the measured quantities and per-check booleans:
        ``p_exc`` [dimensionless], ``p_ground_branch``, ``n_nd``,
        ``raising_max_rel``, ``diag_max_rel``, ``lowering_weight_frac``,
        ``roundtrip_err_rel``, ``unitarity_err``,
        ``relaxation_ok``, ``lowering_only_ok``, ``roundtrip_ok``, ``passed``.
    """
    H, _ = wrap_cqed_hamiltonian(
        qpd, offset_charge, coupling_g_hz, readout_freq_hz, parity=parity,
        n_qubit=n_qubit, n_fock=n_fock, charge_cutoff=charge_cutoff)
    a = cavity_operator(n_qubit, n_fock)
    ops = dressed_jump_operators(
        qpd, offset_charge, coupling_g_hz, readout_freq_hz, parity=parity,
        n_qubit=n_qubit, n_fock=n_fock, charge_cutoff=charge_cutoff)
    U = np.asarray(ops["evecs"])                  # columns = dressed eigvecs
    labels = ops["labels"]                        # (qubit i, photon n) per state
    J_dressed = np.asarray(ops["J_minus_dressed"].full())
    J_bare = np.asarray(ops["J_minus_bare"].full())

    # --- (3) ROUND-TRIP: bare -> dressed -> bare bookkeeping ----------------
    unitarity_err = float(np.max(np.abs(U.conj().T @ U - np.eye(U.shape[0]))))
    Jd_full = U.conj().T @ J_bare @ U             # bare op in dressed basis
    J_rt = U @ Jd_full @ U.conj().T               # ... and back
    scale_bare = float(np.max(np.abs(J_bare)))    # = 1 (unit shift)
    roundtrip_err_rel = float(np.max(np.abs(J_rt - J_bare)) / scale_bare)
    roundtrip_ok = bool(roundtrip_err_rel < roundtrip_tol
                        and unitarity_err < roundtrip_tol)

    # --- (2) ENERGY-LOWERING ONLY (E_k < E_l mask; evals ascending) ---------
    Jd = U.conj().T @ J_dressed @ U               # returned op, dressed basis
    scale_d = float(np.max(np.abs(Jd)))
    tril_strict = np.tril(Jd, k=-1)               # energy-RAISING entries
    diag_part = np.diag(np.diag(Jd))              # energy-conserving entries
    raising_max_rel = float(np.max(np.abs(tril_strict)) / scale_d)
    diag_max_rel = float(np.max(np.abs(diag_part)) / scale_d)
    # Near-null detector (the OLD k=-1 bug signature): the energy-lowering
    # mask must retain essentially all of the projected operator's weight.
    lowering_weight_frac = float(
        np.linalg.norm(Jd, "fro") / np.linalg.norm(Jd_full, "fro"))
    lowering_only_ok = bool(raising_max_rel < mask_noise_tol
                            and diag_max_rel < mask_noise_tol
                            and lowering_weight_frac >= lowering_weight_min)

    # --- (1) RELAXATION: no-drive steady state, J_-^dressed the only qubit
    #         relaxation path (kappa D[a] + Gamma_r D[J_-^dressed]) ----------
    c_ops = [collapse(kappa_hz, a),
             collapse(gamma_r_hz, ops["J_minus_dressed"])]
    rho_ss = steadystate(
        H.to("CSR"), [c.to("CSR") for c in c_ops],
        method="direct", solver="spsolve")
    rho_arr = np.asarray(rho_ss.full())
    pops = np.real(np.diag(U.conj().T @ rho_arr @ U))  # dressed populations
    excited = np.array([lab[0] != 0 for lab in labels], dtype=bool)
    p_exc = float(np.sum(pops[excited]))               # [dimensionless]
    p_ground_branch = float(np.sum(pops[~excited]))
    n_nd = float(np.real(expect(a.dag() * a, rho_ss)))  # cross-ref vs 5.28e-6
    relaxation_ok = bool(p_exc < pexc_tol)

    passed = bool(relaxation_ok and lowering_only_ok and roundtrip_ok)
    result = {
        "n_g": float(offset_charge), "parity": str(parity),
        "n_qubit": int(n_qubit), "n_fock": int(n_fock),
        "p_exc": p_exc, "p_ground_branch": p_ground_branch, "n_nd": n_nd,
        "pexc_tol": float(pexc_tol),
        "raising_max_rel": raising_max_rel, "diag_max_rel": diag_max_rel,
        "lowering_weight_frac": lowering_weight_frac,
        "mask_noise_tol": float(mask_noise_tol),
        "lowering_weight_min": float(lowering_weight_min),
        "roundtrip_err_rel": roundtrip_err_rel,
        "unitarity_err": unitarity_err,
        "roundtrip_tol": float(roundtrip_tol),
        "relaxation_ok": relaxation_ok,
        "lowering_only_ok": lowering_only_ok,
        "roundtrip_ok": roundtrip_ok,
        "passed": passed,
    }
    if assert_pass and not passed:
        raise ValueError(
            "Dressed J_- GKSL round-trip revalidation FAILED at "
            f"n_g={offset_charge} {parity}: relaxation_ok={relaxation_ok} "
            f"(P_exc={p_exc:.3e} vs {pexc_tol:.0e}), "
            f"lowering_only_ok={lowering_only_ok} "
            f"(raising {raising_max_rel:.2e}, diag {diag_max_rel:.2e}, "
            f"weight {lowering_weight_frac:.6f}), "
            f"roundtrip_ok={roundtrip_ok} ({roundtrip_err_rel:.2e}, "
            f"unitarity {unitarity_err:.2e}).")
    return result


# ===========================================================================
# Phase 1.1 Plan 01: LAB-FRAME Floquet driven-dissipative steady state
# ===========================================================================
#
# Assemble the joint H wrap and C5 dressed dissipators with the LAB-FRAME
# cavity drive, solve the PERIODIC (Floquet) steady state once per probe
# frequency omega_d via the ``steadystate_fourier`` Fourier-mode recursion,
# and read S21(omega_d) off the n=+1 sideband field <a>(omega_d) through the
# C2 input-output map.  The decisive observable is the S21 PEAK POSITION
# f_peak (= omega_r + observed-chi), the kappa-insensitive, magnitude-free
# quantity VALD-01 compares to PR #19's compute_observed_chi_spectrum.
#
# LAB FRAME (Phase-1.1 model-form lock; supersedes the Phase-1 rotating
# frame).  In the lab frame the FULL non-RWA coupling g(n_hat-n_g)(a+a_dag)
# is STATIC AND EXACT; the only time dependence is the monochromatic cavity
# drive:
#
#     H_lab(t) = H_cQED + 2*eps*cos(omega_d t)*(a + a_dag)
#     H_cQED   = diag(omega_q_i) (x) I_r + omega_r I_q (x) a_dag a
#                + g(n_hat - n_g)(a + a_dag)          [static, non-RWA, exact]
#
# There is NO -omega_d*a_dag a term anywhere (fp-rf-shortcut guard): holding
# the non-RWA coupling static in the omega_d rotating frame was the indicted
# Phase-1 VALD-01 bug (g*n_hat*a_dag became a spurious static cavity source,
# no-drive <n>~1.53, sign-flipped peak).
#
# DRIVE NORMALIZATION (C3 <-> lab frame).  The convention lock defines eps in
# the ROTATING-FRAME form H_d^RF = eps*(a+a_dag) with alpha_in = eps/sqrt(
# kappa_c) and the SIMU-01 closed form <a>(omega_r) = -2i*eps/kappa.  The
# lab-frame drive that realizes the SAME eps is
#
#     H_d(t) = 2*eps*cos(omega_d t)*(a+a_dag)
#            = eps*(e^{+i omega_d t} + e^{-i omega_d t})*(a+a_dag),
#
# i.e. each e^{-+i omega_d t} Fourier sideband carries the full eps*(a+a_dag)
# (its RWA/rotating-frame limit is exactly the C3-locked eps*(a+a_dag)).
# With the half-amplitude choice (eps*cos) the g=0 sideband comes out
# -i*eps/kappa -- off the SIMU-01 closed form by exactly 2 and the notch
# S21(omega_r) = 1-2*kappa_c/kappa breaks (verified numerically, Plan-01.1-01
# DEVIATION Rule 4) -- so the factor of 2 belongs to the drive operator, NOT
# to s21/alpha_in, which are reused verbatim.
#
# FOURIER-MODE (Floquet) STEADY STATE.  The periodic state is
# rho(t) = sum_n rho_n e^{-i n omega_d t}.  QuTiP's ``steadystate_fourier``
# (5.3; ``steadystate_floquet`` is the deprecated alias) solves the coupled
# harmonic system with the Sze Meng Tan matrix-continued-fraction recursion
# (quantum optics toolbox section 16) and returns the DC component rho_0.
# It makes NO secular/Markov approximation on the system: the full non-RWA
# coupling and the dressed dissipators enter L_0 = liouvillian(H_cQED, c_ops)
# verbatim; only the periodic structure is expanded (truncated at n_it
# harmonics).  Floquet-Markov machinery (fmmesolve/FloquetBasis/
# floquet_tensor) is FORBIDDEN here (fp-fmmesolve): its secular step is a
# hidden RWA-like reduction.
#
# DC-vs-SIDEBAND (the central trap, fp-dc-sideband).  rho_0 is the TIME-
# AVERAGE: Tr(a rho_0) = 0 exactly (the field oscillates at omega_d).  The
# readout field lives in the n=+1 sideband
#
#     rho_+1 = vector_to_operator(T @ vec(rho_0)),
#     <a>(omega_d) = Tr(a rho_+1)        [the e^{-i omega_d t} component]
#
# while the DC occupation <n> = Tr(a_dag a rho_0) is the true weak-drive
# cavity occupation.  The public steadystate_fourier returns ONLY rho_0 (no
# sideband accessor -- confirmed by source inspection), so the recursion is
# reimplemented below, CSR end-to-end, ALSO returning the S/T sideband maps.
# Sideband identity VERIFIED on the g=0 cavity: Tr(a rho_+1) reproduces the
# rotating-frame Lorentzian -i*eps_sb/(i*(omega_r-omega_d)+kappa/2) exactly
# (eps_sb the per-sideband amplitude), while Tr(a rho_-1) ~ 0 (rho_-1 carries
# <a_dag>).
#
# The collapse set is the FULL dressed set (C5), unchanged:
#     [ collapse(kappa, a),               # cavity, bare field, TOTAL kappa
#       collapse(Gamma_r, J_-^dressed),   # dressed qubit relaxation
#       collapse(Gamma_phi/2, Sz_dressed) ]  # dressed pure dephasing
#
# STEADY-STATE UNIQUENESS (carried over from Phase 1, still true in the lab
# frame).  A qubit-relaxation channel Gamma_r > 0 is REQUIRED: with the
# cavity dissipator alone the joint Liouvillian is singular (the qubit has no
# relaxation path) and the linear solves do not converge.  Gamma_r only
# selects the |0,n> ground branch and adds sub-kappa broadening; the peak
# position is Gamma_r-insensitive (fp-magnitude honored: kappa and rates are
# never tuning levers).
# ===========================================================================


# ===========================================================================
# Phase 2 Plan 02-01: DISPLACED-FRAME solver mode (DERV-02)
# ===========================================================================
#
# Lab-frame displacement a -> a + alpha(t) with the bare-cavity classical
# alpha(t) solving (R.1)  alpha_dot = -(i*omega_r + kappa/2)*alpha
# - i*2*eps*cos(omega_d t)  (DERV-02 section 2; the kappa/2 is DERIVED from the
# kappa*D[a] linear drift).  The unique periodic solution is (R.3)
#
#     alpha(t) = alpha_+ e^{-i omega_d t} + alpha_- e^{+i omega_d t}
#     alpha_+  = -i*eps / [kappa/2 + i*(omega_r - omega_d)]
#     alpha_-  = -i*eps / [kappa/2 + i*(omega_r + omega_d)]   (counter-rotating,
#                                                              kept exactly)
#
# In the displaced (residual) frame the cavity drive is GONE and the drive
# moves to the qubit charge operator (R.2):
#
#     H_res(t) = H_cQED + g(n_hat - n_g) * 2*Re alpha(t)
#              = H_cQED + V_+ e^{-i omega_d t} + V_- e^{+i omega_d t}
#     V_+ = g(n_hat - n_g) * A_+,   A_+ = alpha_+ + conj(alpha_-)
#     V_- = g(n_hat - n_g) * A_-,   A_- = conj(A_+)   =>   V_- = V_+^dag
#
# L0 (H_cQED + kappa*D[a] + dressed J_-/Sz, C5) is UNCHANGED, so the same
# truncated Fourier recursion applies with UNEQUAL sideband weights.  The
# legacy solver is the special case V_+ = V_- = eps*(a+a_dag) at alpha = 0
# (the bit-for-bit backward-compatibility default).
#
# HARMONIC SLOT MAPPING (load-bearing for unequal weights).  With
# rho(t) = sum_n rho_n e^{-i n omega_d t} and
# L(t) = L0 + L_plus e^{-i omega_d t} + L_minus e^{+i omega_d t}, matching
# e^{-i n omega_d t} components gives
#
#     T_n = -(L0 + i n w + L_minus T_{n+1})^{-1} L_plus      (rho_n = T_n rho_{n-1})
#     S_n = -(L0 - i n w + L_plus  S_{n+1})^{-1} L_minus     (rho_{-n} = S_n rho_{-(n-1)})
#
# i.e. the legacy variable named ``Lm`` is the L_plus slot and the legacy
# ``Lp`` is the L_minus slot (historically swapped names; invisible at equal
# weights, pinned here by the g=0 anchor -2i*eps/kappa and the fast-vs-
# production cross-validation in checks/displaced_anchor_checks.py).
#
# FAST-PATH NEGATIVE-CHAIN SYMMETRY (DERV-02 Appendix A).  The reuse
# C_minus = conj(C_plus[perm, perm]) and rho_{-n} = rho_n^dag holds iff
# V_- = V_+^dag (Hermitian H_res(t)); it does NOT require real or equal
# weights.  The displaced mode satisfies V_- = V_+^dag by construction and
# ASSERTS it (plus the matrix-level P conj(L_G) P = L_G identity) at
# construction time; on failure the fast path refuses to run (the production
# solver, which computes the S chain independently, is the full-chain
# fallback).
#
# READOUT RECONSTRUCTION (exact, DERV-02 section 5; C2 map unchanged):
#
#     <a>_lab(omega_d) = alpha_+ + Tr(a rho_+1^res)          (5.1)
#     <n>_lab = Tr(a_dag a rho_0^res)
#             + 2*Re[conj(alpha_+) Tr(a rho_+1^res)]
#             + 2*Re[conj(alpha_-) Tr(a rho_-1^res)]
#             + |alpha_+|^2 + |alpha_-|^2                     (5.3)
#
# Leakage observables (Goto-Koshino-style monitors, arXiv:2305.00628):
# p_top (top-Fock population of the residual state, max over stored
# sidebands), <n>_res, |Tr(a rho_+1^res)|.
# ===========================================================================


def displacement_alphas(f_r_hz: float, omega_d_hz: float, eps_hz: float,
                        kappa_hz: float):
    """Closed-form alpha_+/- of the periodic bare-cavity displacement (R.3).

    ALL-Hz RATIO (single-2*pi-site discipline, C1): alpha_+- are dimensionless
    ratios of frequency-dimensioned quantities, so the 2*pi cancels and the
    formula is evaluated directly on Hz API quantities -- no second conversion
    site is created.

        alpha_+ = -i*eps_hz / [kappa_hz/2 + i*(f_r - f_d)]
        alpha_- = -i*eps_hz / [kappa_hz/2 + i*(f_r + f_d)]

    Signs/phases are pinned against the C2 -i drive phase by the g=0 anchor:
    alpha_+(omega_d = omega_r) = -2i*eps/kappa (the locked DERV-01 sideband
    value, rel err 8.3e-13).  The counter-rotating alpha_- is kept exactly
    (non-RWA policy; |alpha_-| ~ eps/(2*omega_r)).

    Returns
    -------
    (alpha_p, alpha_m) : complex, complex
        The e^{-i omega_d t} and e^{+i omega_d t} coefficients [dimensionless].
    """
    alpha_p = -1j * eps_hz / (0.5 * kappa_hz + 1j * (f_r_hz - omega_d_hz))
    alpha_m = -1j * eps_hz / (0.5 * kappa_hz + 1j * (f_r_hz + omega_d_hz))
    return complex(alpha_p), complex(alpha_m)


def charge_drive_operator(components: dict, coupling_g_hz: float) -> Qobj:
    """The residual-drive operator g*(n_hat - n_g) (x) I_r in Hz (DERV-02 R.2).

    Built from the ``n_shift_mat`` component returned by
    :func:`wrap_cqed_hamiltonian` (the n_hat - n_g_eff matrix in the truncated
    qubit eigenbasis, parity shift already applied), qubit-major (C1).  The
    returned operator is in Hz; the single 2*pi site (``hz_to_qobj``) is
    applied INSIDE the displaced solvers, never here.
    """
    n_shift = np.asarray(components["n_shift_mat"], dtype=float)
    n_qubit = int(components["n_qubit"])
    n_fock = int(components["n_fock"])
    dims = [[n_qubit, n_fock], [n_qubit, n_fock]]
    return Qobj(coupling_g_hz * np.kron(n_shift, np.eye(n_fock)), dims=dims)


def _displaced_sideband_weights(f_r_hz, omega_d_hz, eps_hz, kappa_hz):
    """(alpha_p, alpha_m, A_p) with A_p = alpha_p + conj(alpha_m) (DERV-02 4.1)."""
    alpha_p, alpha_m = displacement_alphas(f_r_hz, omega_d_hz, eps_hz, kappa_hz)
    return alpha_p, alpha_m, alpha_p + np.conj(alpha_m)


def _check_drive_op_hermitian(drive_op_hz: Qobj, tol: float = 1e-12):
    """Construction-time guard: V_- = V_+^dag requires a Hermitian charge op.

    DERV-02 Appendix A failure mode: a non-Hermitian drive operator breaks
    the negative-chain symmetry SILENTLY -- assert, never assume.
    """
    herm_rel = (drive_op_hz - drive_op_hz.dag()).norm() / max(
        drive_op_hz.norm(), 1.0)
    if herm_rel > tol:
        raise ValueError(
            "Displaced-mode drive operator is not Hermitian "
            f"(rel err {herm_rel:.2e} > {tol:.0e}); V_- = V_+^dag would fail "
            "and the fast-path negative-chain symmetry (DERV-02 App. A) "
            "would silently solve the wrong system.")


def _top_fock_trace(rho_arr: np.ndarray, n_qubit: int, n_fock: int) -> complex:
    """sum_q <q, N_f-1| rho |q, N_f-1> for a qubit-major joint matrix (C1)."""
    idx = np.arange(n_qubit) * n_fock + (n_fock - 1)
    return complex(np.sum(rho_arr[idx, idx]))


def driven_steadystate(
    H_cqed,
    a,
    omega_d_hz: float,
    eps_hz: float,
    kappa_hz: float,
    c_ops_extra=None,
    n_it: int = 3,
    displaced: bool = False,
    f_r_hz: float = None,
    drive_op_hz: Qobj = None,
):
    """Lab-frame Floquet (periodic) steady state at probe frequency omega_d.

    Solves the periodic steady state of

        H_lab(t) = H_cQED + 2*eps*cos(omega_d t)*(a + a_dag)

    with the ``steadystate_fourier`` Fourier-mode recursion (CSR/spsolve,
    reimplemented so the n=+/-1 sideband maps S, T are available) and returns
    BOTH the readout sideband field ``a_wd`` = <a>(omega_d) = Tr(a rho_+1)
    and the DC occupation ``n_dc`` = Tr(a_dag a rho_0).  The two live in
    DIFFERENT Fourier components -- rho_0 is the time-average (its <a> is 0
    by construction), rho_+1 carries the field at omega_d (fp-dc-sideband).

    ``H_cqed`` enters STATIC and UNSHIFTED: there is NO -omega_d*a_dag a term
    (fp-rf-shortcut -- the indicted Phase-1 rotating-frame bug).  The factor
    2 on the drive makes each e^{-+i omega_d t} sideband carry the C3-locked
    eps*(a+a_dag), so the SIMU-01 closed forms (<a>(omega_r) = -2i*eps/kappa
    at g=0) and the verbatim ``s21`` map hold unchanged (see section note).

    All frequency inputs are in Hz; the SINGLE 2*pi site is ``hz_to_qobj``
    (drive operator and omega_d alike) and the collapse rate carries 2*pi via
    ``collapse`` (C1).

    Parameters
    ----------
    H_cqed : Qobj
        Joint qubit+cavity Hamiltonian in rad/s (from ``wrap_cqed_hamiltonian``);
        carries ``omega_r a_dag a``, the multilevel transmon, and the non-RWA
        coupling ``g(n_hat-n_g)(a+a_dag)``.  Used AS IS (lab frame).
    a : Qobj
        Bare cavity annihilation in the joint space (``cavity_operator``).
    omega_d_hz : float
        Probe/drive frequency omega_d [Hz].
    eps_hz : float
        Drive amplitude epsilon [Hz] in the C3 rotating-frame normalization
        (alpha_in = eps/sqrt(kappa_c)); the lab-frame cos drive is 2*eps.
    kappa_hz : float
        TOTAL cavity linewidth kappa [Hz] (the cavity dissipator rate, C2/C5).
    c_ops_extra : list of Qobj, optional
        Extra (already 2*pi-scaled) collapse operators, e.g.
        ``collapse(Gamma_r_hz, J_minus_dressed)`` and
        ``collapse(Gamma_phi_hz/2, Sz_dressed)``.  At least one dressed qubit
        relaxation term (Gamma_r > 0) is REQUIRED for a unique steady state
        (see section note).
    n_it : int
        Number of +/- Fourier harmonics kept in the recursion.  n_it=3
        resolves the +/- sideband asymmetry at weak drive (<n> ~ 0.01); bump
        3->5 as the Plan-02 convergence check.  Grows with drive power
        (Phase-3 concern).
    displaced : bool
        If True, solve the DISPLACED-FRAME residual problem (DERV-02 R.2):
        the cavity drive is removed exactly by the classical bare-cavity
        alpha(t) of (R.3) and the sideband operators become
        V_+- = g(n_hat-n_g)*A_+- (unequal weights, V_- = V_+^dag).  The
        returned ``a_wd``/``n_dc`` are then the EXACT lab-frame
        reconstructions (5.1)/(5.3) -- the C2 ``s21`` map applies verbatim.
        Default False reproduces the legacy solver bit-for-bit.
    f_r_hz : float, required if displaced
        Bare resonator frequency omega_r [Hz] entering (R.3).
    drive_op_hz : Qobj, required if displaced
        The Hz-dimensioned residual drive operator g*(n_hat-n_g) (x) I_r in
        the joint space (``charge_drive_operator``); must be Hermitian
        (asserted -- DERV-02 App. A guard).

    Returns
    -------
    dict with keys
        ``a_wd``   : complex, the readout field <a>(omega_d) [dimensionless]
                     -- feed THIS into ``s21``.  Default mode: Tr(a rho_+1).
                     Displaced mode: alpha_+ + Tr(a rho_+1^res) (exact lab
                     reconstruction, DERV-02 (5.1));
        ``n_dc``   : float, cavity occupation [dimensionless].  Default mode:
                     Tr(a_dag a rho_0).  Displaced mode: the exact <n>_lab
                     reconstruction (5.3) incl. cross terms and |alpha_+-|^2;
        ``a_dc``   : complex, DC field Tr(a rho_0) (~0 by construction; kept
                     as the fp-dc-sideband self-check);
        ``rho_0``  : Qobj, the DC (time-averaged) density matrix (residual-
                     frame in displaced mode);
        ``rho_p1`` : Qobj, the n=+1 sideband component (NOT a density matrix:
                     traceless, non-Hermitian; it is a Fourier coefficient);
        ``rho_m1`` : Qobj, the n=-1 sideband (= rho_p1^dag by the App. A
                     Hermitian-pair identity; computed INDEPENDENTLY from the
                     S chain here, so comparing it against rho_p1.dag() is an
                     executable check of that identity);
        ``a_res``  : complex, residual coherent amplitude Tr(a rho_+1)
                     (obs-res-amp; equals a_wd in default mode);
        ``n_res``  : float, residual occupation Tr(a_dag a rho_0)
                     (equals n_dc in default mode);
        ``p_top``  : float, top-Fock population of rho_0 (truncation monitor);
        ``p_top_max`` : float, max |top-Fock trace| over stored sidebands;
        ``alpha_p``, ``alpha_m`` : complex, the (R.3) displacement
                     coefficients (0 in default mode).
    """
    # --- Lab-frame pieces (NO omega_d shift on H; fp-rf-shortcut guard) -----
    c_ops = [collapse(kappa_hz, a)]
    if c_ops_extra:
        c_ops = c_ops + list(c_ops_extra)

    w_d = hz_to_qobj(omega_d_hz)  # 2*pi*omega_d [rad/s], single 2*pi site

    # --- steadystate_fourier recursion (QuTiP 5.3 semantics), CSR/spsolve ---
    # rho(t) = sum_n rho_n e^{-i n omega_d t}; downward matrix-continued-
    # fraction recursion over n_it harmonics.  Everything is kept CSR and
    # every linear solve goes through the sparse LU path ('spsolve') -- the
    # joint Liouvillian is (n_q*n_f)^2 a side and dense solves are
    # prohibitive (cf. the Plan-02 Purcell note).
    #
    # SLOT MAPPING (section note above): legacy ``Lm`` carries L_plus (the
    # e^{-i omega_d t} sideband generator) and legacy ``Lp`` carries L_minus.
    L0 = liouvillian(H_cqed, c_ops).to("CSR")
    if displaced:
        if f_r_hz is None or drive_op_hz is None:
            raise ValueError(
                "displaced=True requires f_r_hz and drive_op_hz "
                "(charge_drive_operator).")
        _check_drive_op_hermitian(drive_op_hz)
        alpha_p, alpha_m, A_p = _displaced_sideband_weights(
            f_r_hz, omega_d_hz, eps_hz, kappa_hz)
        G_rad = hz_to_qobj(drive_op_hz)        # single 2*pi site (C1)
        V_plus = A_p * G_rad
        V_minus = np.conj(A_p) * G_rad
        # Construction-time assertion of the App. A condition V_- = V_+^dag.
        vnorm = V_plus.norm()
        if vnorm > 0.0 and (V_minus - V_plus.dag()).norm() / vnorm > 1e-12:
            raise ValueError(
                "Displaced sideband operators violate V_- = V_+^dag; "
                "negative-chain symmetry (DERV-02 App. A) does not hold.")
        # CRITICAL: the sideband generator is the COMMUTATOR superoperator
        # -i[V_pm, .] = A_pm * (-i[G, .]).  QuTiP's liouvillian(V) for a
        # NON-Hermitian V builds the effective-Hamiltonian generator
        # -i(V rho - rho V^dag) instead -- which is trace-non-conserving and
        # silently breaks the harmonic hierarchy (caught executably: the DC
        # operator loses its left trace-null vector and rho_-1 != rho_+1^dag).
        # Therefore scale the HERMITIAN-G commutator Liouvillian by the
        # complex weights; NEVER call liouvillian(A_pm * G_rad) directly.
        LG = liouvillian(G_rad).to("CSR")      # -i[G, .]  (G Hermitian)
        Lm = (A_p * LG).to("CSR")              # L_plus  (e^{-i w_d t})
        Lp = (np.conj(A_p) * LG).to("CSR")     # L_minus (e^{+i w_d t})
    else:
        # Legacy equal-weight cavity drive: H_d(t) = Op_t*cos(omega_d t),
        # Op_t = 2*eps*(a+a_dag) so each sideband is the C3 eps*(a+a_dag).
        # This branch is BYTE-IDENTICAL to the pre-displacement solver (the
        # alpha=0 backward-compatibility toggle, checks/ anchor (b)).
        alpha_p = 0.0 + 0.0j
        alpha_m = 0.0 + 0.0j
        Op_t = hz_to_qobj(2.0 * eps_hz) * (a + a.dag())
        Lm = (0.5 * liouvillian(Op_t)).to("CSR")
        Lp = (0.5 * liouvillian(Op_t)).to("CSR")
    Id = qeye_like(L0).to("CSR")
    S = qzero_like(L0).to("CSR")
    T = qzero_like(L0).to("CSR")
    for n_i in range(n_it, 0, -1):
        L = (L0 - 1j * n_i * w_d * Id + Lm @ S).to("CSR")
        S.data = -_data.solve(L.data, Lp.data, "spsolve")
        S = S.to("CSR")
        L = (L0 + 1j * n_i * w_d * Id + Lp @ T).to("CSR")
        T.data = -_data.solve(L.data, Lm.data, "spsolve")
        T = T.to("CSR")

    # DC component rho_0 from the harmonic-dressed Liouvillian.
    rho_0 = steadystate((L0 + Lm @ S + Lp @ T).to("CSR"), solver="spsolve")

    # n=+1 sideband (the e^{-i omega_d t} component): carries <a>(omega_d).
    rho_p1 = vector_to_operator(T @ operator_to_vector(rho_0))
    # n=-1 sideband from the INDEPENDENT S chain (carries <a_dag>); must equal
    # rho_p1^dag by the Hermitian-pair identity (DERV-02 App. A) -- returned
    # so the anchor checks can verify that identity executably.
    rho_m1 = vector_to_operator(S @ operator_to_vector(rho_0))

    a_res = complex(expect(a, rho_p1))                 # residual/readout field
    n_res = float(np.real(expect(a.dag() * a, rho_0)))  # residual occupation
    a_dc = complex(expect(a, rho_0))                   # ~0 (time-avg) self-check

    # Leakage observables (DERV-02 section 5.3): top-Fock truncation monitor.
    n_qubit, n_fock = a.dims[0]
    p_top = float(np.real(_top_fock_trace(rho_0.full(), n_qubit, n_fock)))
    p_top_max = max(
        abs(_top_fock_trace(rho_0.full(), n_qubit, n_fock)),
        abs(_top_fock_trace(rho_p1.full(), n_qubit, n_fock)),
        abs(_top_fock_trace(rho_m1.full(), n_qubit, n_fock)),
    )

    if displaced:
        # Exact lab-frame readout reconstruction (DERV-02 (5.1), (5.3)).
        a_wd = alpha_p + a_res
        a_res_m1 = complex(expect(a, rho_m1))          # Tr(a rho_-1^res)
        n_dc = (n_res
                + 2.0 * float(np.real(np.conj(alpha_p) * a_res))
                + 2.0 * float(np.real(np.conj(alpha_m) * a_res_m1))
                + abs(alpha_p) ** 2 + abs(alpha_m) ** 2)
    else:
        a_wd = a_res
        n_dc = n_res

    return {
        "a_wd": a_wd,
        "n_dc": n_dc,
        "a_dc": a_dc,
        "rho_0": rho_0,
        "rho_p1": rho_p1,
        "rho_m1": rho_m1,
        "a_res": a_res,
        "n_res": n_res,
        "p_top": p_top,
        "p_top_max": float(p_top_max),
        "alpha_p": alpha_p,
        "alpha_m": alpha_m,
    }


def nodrive_precondition(
    H_cqed,
    a,
    kappa_hz: float,
    c_ops_extra=None,
    tol_floor: float = 1e-3,
    assert_floor: bool = True,
):
    """No-drive <n> -> 0 hard precondition (eps = 0 vacuum check).

    Solves the UNDRIVEN steady state of the static lab-frame Hamiltonian
    ``H_lab = H_cQED`` (NO drive, NO -omega_d*a_dag a shift -- subtracting one
    is fp-rf-shortcut, the indicted Phase-1 rotating-frame bug) with the full
    C5 dressed collapse set, and checks that the cavity relaxes to (dressed)
    vacuum:

        rho_nd = steadystate(H_lab, c_ops)  at eps = 0
        <n>_no-drive = Tr(a_dag a rho_nd)  <  tol_floor

    Physical bound.  The dressed ground state carries a bare-photon admixture
    of order (g_01/Delta)^2 ~ (27.9 MHz / 1.622 GHz)^2 ~ 3e-4 at the QPD4
    point (g_01 = g*|n_01|, n_01 = 0.6968; Delta = -1.622 GHz), so
    ``tol_floor = 1e-3`` is a generous bound on the physical floor.  The
    verified value at QPD4 n_g=0 even (n_qubit=8, n_fock=6,
    Gamma_r = kappa/100) is ~5.6e-6 -- six orders of magnitude below the
    unphysical no-drive <n> ~ 1.53 of the rotating-frame static-non-RWA bug
    (the VALD-01 HALT root cause).  A no-drive <n> of O(1) here is the single
    cleanest disconfirmation that fp-rf-shortcut has regressed in.

    All rates are in Hz; the collapse rate carries the 2*pi via ``collapse``
    (C1).  The solve is CSR/spsolve end-to-end (RESEARCH.md "No-drive
    precondition wiring").

    Parameters
    ----------
    H_cqed : Qobj
        Static joint Hamiltonian in rad/s (``wrap_cqed_hamiltonian``), used AS
        IS -- lab frame, no omega_d shift.
    a : Qobj
        Bare cavity annihilation in the joint space (``cavity_operator``).
    kappa_hz : float
        TOTAL cavity linewidth kappa [Hz] (cavity dissipator rate, C2/C5).
    c_ops_extra : list of Qobj, optional
        Extra (already 2*pi-scaled) collapse operators -- the dressed qubit
        relaxation/dephasing set.  Gamma_r > 0 is REQUIRED for a unique
        steady state (same uniqueness note as :func:`driven_steadystate`).
    tol_floor : float
        Vacuum-floor bound on <a_dag a> (default 1e-3, generous vs the
        ~3e-4 dressed-admixture scale).
    assert_floor : bool
        If True (default), raise ``ValueError`` when the floor is violated;
        if False, only report (the caller gates).

    Returns
    -------
    dict with keys
        ``n_nd`` : float, the no-drive occupation <a_dag a> [dimensionless];
        ``passed`` : bool, ``n_nd < tol_floor``;
        ``tol_floor`` : float, the bound used;
        ``rho_nd`` : Qobj, the no-drive steady state.
    """
    c_ops = [collapse(kappa_hz, a)]
    if c_ops_extra:
        c_ops = c_ops + list(c_ops_extra)
    rho_nd = steadystate(
        H_cqed.to("CSR"),
        [c.to("CSR") for c in c_ops],
        method="direct",
        solver="spsolve",
    )
    n_nd = float(np.real(expect(a.dag() * a, rho_nd)))
    passed = n_nd < tol_floor
    if assert_floor and not passed:
        raise ValueError(
            f"No-drive precondition FAILED: <a_dag a>_no-drive = {n_nd:.3e} "
            f">= tol_floor = {tol_floor:.1e}. An O(1) value signals the "
            "rotating-frame static-non-RWA shortcut (fp-rf-shortcut) has "
            "regressed in; expected ~5.6e-6 at the QPD4 point."
        )
    return {"n_nd": n_nd, "passed": passed, "tol_floor": tol_floor, "rho_nd": rho_nd}


def driven_s21_sweep(
    H_cqed,
    a,
    omega_d_grid_hz,
    eps_hz: float,
    kappa_hz: float,
    kappa_c_hz: float,
    c_ops_extra=None,
    n_it: int = 3,
):
    """Sweep omega_d to build S21(omega_d) from the lab-frame Floquet sideband.

    Calls :func:`driven_steadystate` (lab-frame Floquet) once per omega_d in
    ``omega_d_grid_hz`` and maps the SIDEBAND field ``a_wd`` = Tr(a rho_+1)
    -- NOT the DC <a>, which is ~0 by construction (fp-dc-sideband) --
    through the C2 input-output relation
    ``S21 = 1 - i*sqrt(kappa_c)<a>/alpha_in`` with ``alpha_in = eps/sqrt(kappa_c)``
    (so ``S21 = 1 - i*(kappa_c/eps)<a>`` -- note <a> is itself linear in eps at
    weak drive, so S21 is drive-amplitude independent in the linear-response
    corner; this is the verify check).

    Parameters
    ----------
    omega_d_grid_hz : array_like
        Probe frequencies [Hz] (typically a window +/- W*kappa about the
        dressed cavity).
    kappa_hz : float
        TOTAL kappa [Hz] -- the dissipator rate.
    kappa_c_hz : float
        PORT kappa_c [Hz] -- the S21 prefactor rate (C2; kept distinct from
        total kappa).
    n_it : int
        Floquet harmonic truncation passed to :func:`driven_steadystate`.

    Returns
    -------
    dict with keys
        ``omega_d_hz`` (the grid), ``a_expect`` (complex SIDEBAND field
        <a>(omega_d) = Tr(a rho_+1)), ``s21`` (complex S21(omega_d)), and
        ``n_expect`` (the DC cavity occupation <a_dag a>(rho_0), for
        confirming the weak-drive corner <n> << 1).
    """
    grid = np.asarray(omega_d_grid_hz, dtype=float)
    a_vals = np.empty(grid.size, dtype=complex)
    n_vals = np.empty(grid.size, dtype=float)
    for j, wd in enumerate(grid):
        res = driven_steadystate(
            H_cqed, a, wd, eps_hz, kappa_hz,
            c_ops_extra=c_ops_extra, n_it=n_it,
        )
        a_vals[j] = res["a_wd"]   # sideband field (rho_+1), NOT DC <a>
        n_vals[j] = res["n_dc"]   # DC occupation (rho_0)
    s21_vals = np.array([s21(av, eps_hz, kappa_c_hz) for av in a_vals])
    return {
        "omega_d_hz": grid,
        "a_expect": a_vals,
        "s21": s21_vals,
        "n_expect": n_vals,
    }


def extract_peak_hz(omega_d_grid_hz, s21_vals):
    """Sub-linewidth S21 peak position by parabolic interpolation of |S21 - 1|.

    The transmission notch/peak is where the cavity response |S21 - 1| (the
    deviation from full transmission) is LARGEST -- i.e. the intracavity field
    <a> is resonantly enhanced.  We locate the discrete argmax of |S21 - 1| then
    refine to sub-grid resolution by a 3-point parabolic fit (vertex of the
    parabola through the peak sample and its two neighbours).  This gives a
    peak position accurate to well below one grid step (the resolution that
    matters for the ~0.05*kappa VALD-01 tolerance).

    Parameters
    ----------
    omega_d_grid_hz : array_like
        Monotonic probe-frequency grid [Hz] (uniform spacing assumed for the
        parabolic refinement).
    s21_vals : array_like of complex
        S21(omega_d) on the grid.

    Returns
    -------
    f_peak_hz : float
        Interpolated peak position [Hz].
    """
    grid = np.asarray(omega_d_grid_hz, dtype=float)
    resp = np.abs(np.asarray(s21_vals) - 1.0)  # cavity response = |S21 - 1|
    k = int(np.argmax(resp))
    # Edge guard: if the peak is at a grid edge the window is mis-centred;
    # return the discrete argmax (the caller widens/recentres the window).
    if k == 0 or k == grid.size - 1:
        return float(grid[k])
    # 3-point parabolic vertex (uniform grid): offset in units of the step.
    y0, y1, y2 = resp[k - 1], resp[k], resp[k + 1]
    denom = (y0 - 2.0 * y1 + y2)
    if denom == 0.0:
        return float(grid[k])
    delta = 0.5 * (y0 - y2) / denom  # vertex offset in [-1, 1] grid steps
    step = grid[k + 1] - grid[k]
    return float(grid[k] + delta * step)


# ===========================================================================
# Phase 1.1 Plan 02: TIGHTENED, FACTORED VALD-01 gate harness
# ===========================================================================
#
# The gate is a self-consistency check of ONE Hamiltonian against itself (the
# dressed |0,0> -> |0,1> cavity pole), FACTORED so a model-form bug and a
# convention/benchmark issue localize to different legs:
#
#   PRECONDITION (hard):  no-drive <a_dag a> < tol_floor at each (n_g,parity)
#                         BEFORE any peak comparison (``nodrive_precondition``).
#   LEG (a):  |f_peak^driven - omega_r - chi_numeric[0]| < tol_a, where
#             chi_numeric = compute_dispersive_matrix(method='numeric',
#             rwa=False)[1] on the IDENTICAL build_cqed_hamiltonian params
#             -- isolates MODEL-FORM/solver bugs.
#   LEG (b):  |chi_numeric[0] - chi_obs[0]| < tol_b with crossing[0]==False
#             and weight_frac[0] >= 0.9 asserted first (PR #19
#             compute_observed_chi_spectrum) -- isolates CONVENTION/BENCHMARK
#             issues.  Never the direct driven-pole-vs-chi_obs shortcut
#             (fp-leg-conflation).
#
# The tolerance is MEASURED, not asserted:
#   tol = max(grid_resolution, truncation_drift[nf+1, nq+1, n_it 3->5],
#             stark_term = <n>*dchi/dn).
#
# DENSIFIED GRID + L0 CACHING (Plan-01 carry-forward: ~60 s/point on the CSR
# map recursion makes a ~0.01*kappa grid prohibitive).  The sweep therefore
# uses ``FloquetSidebandSolver``: a dense-LAPACK solver of the IDENTICAL
# truncated Fourier-harmonic system that ``driven_steadystate`` solves.
# Writing rho(t) = sum_n rho_n e^{-i n omega_d t} and matching e^{-i n
# omega_d t} coefficients of d rho/dt = L(t) rho with
# L(t) = L0 + Lc e^{+i omega_d t} + Lc e^{-i omega_d t}, Lc = 0.5*L[Op_t]:
#
#     (L0 + i n w_d) rho_n + Lc rho_{n-1} + Lc rho_{n+1} = 0,  |n| <= n_it,
#     rho_{+-(n_it+1)} = 0,   Tr(rho_0) = 1,
#
# which is algebraically the SAME system the steadystate_fourier downward
# S/T matrix-continued-fraction recursion eliminates (the T chain is
# rho_n = T_n rho_{n-1} for n>0, the S chain rho_{-n} = S_n rho_{-(n-1)};
# substituting recovers the block equations above, same truncation).  The
# fast path is therefore NOT an approximation of the production solver --
# it is the same linear system solved by block-preconditioned GMRES instead
# of dense map continued fractions, and it is cross-validated against
# ``driven_steadystate`` to solver precision in the Plan-02 check script
# before any gate use.  L0 (omega_d-independent) is built and trace-bordered
# LU-factorized ONCE per (n_g,parity) -- the contract's L0-caching item.
# ===========================================================================


class FloquetSidebandSolver:
    """Cached dense solver of the lab-frame Floquet harmonic system.

    Solves, per probe frequency omega_d, the truncated Fourier-harmonic
    steady-state system of ``driven_steadystate`` (see section note above --
    identical physics, identical truncation, different linear algebra) and
    returns the readout sideband field ``a_wd`` = Tr(a rho_+1), the DC
    occupation ``n_dc`` = Tr(a_dag a rho_0), and the DC field self-check
    ``a_dc`` ~ 0 (fp-dc-sideband guard).

    The omega_d-INDEPENDENT pieces are cached at construction (the Plan-02
    L0-caching requirement):

      * ``L0`` = liouvillian(H_cqed, c_ops) as a dense complex matrix,
      * the unit-eps drive Liouvillian ``Lc_unit`` (Lc = eps_hz * Lc_unit;
        Op_t = 2*eps*(a+a_dag) so each e^{-+i omega_d t} sideband carries the
        C3-locked eps*(a+a_dag), Plan-01 factor-2 lock),
      * the trace-bordered LU of L0 (row ``trace_row`` replaced by the trace
        functional) used for the n=0 block and the no-drive steady state.

    Per omega_d only the shifted diagonal blocks (L0 + i n w_d) are
    LU-factorized.  Solves with (L0 - i n w_d) reuse the SAME factorizations
    through the Hermiticity-preservation symmetry of any GKSL generator,
    L0 = P conj(L0) P with P the vec-transpose permutation (verified
    numerically at construction; GKSL maps satisfy L(X^dag) = (L X)^dag).

    Parameters
    ----------
    H_cqed : Qobj
        Static lab-frame joint Hamiltonian in rad/s (``wrap_cqed_hamiltonian``),
        used AS IS -- no -omega_d*a_dag a shift (fp-rf-shortcut guard).
    a : Qobj
        Bare cavity annihilation in the joint space (``cavity_operator``).
    kappa_hz : float
        TOTAL cavity linewidth [Hz] (cavity dissipator rate, C2/C5).
    c_ops_extra : list of Qobj, optional
        Extra (already 2*pi-scaled) collapse operators -- the C5 dressed set.
    n_it : int
        Default Fourier-harmonic truncation (same meaning as
        ``driven_steadystate``); overridable per solve for the convergence bump.
    displaced : bool
        If True, run the DISPLACED-FRAME residual mode (DERV-02 R.2): per
        solve the sideband generators are L_+- = A_+- * liouvillian(G) with
        G = 2*pi*g(n_hat-n_g) and A_- = conj(A_+) from (R.3)/(4.1); the
        returned ``a_wd``/``n_dc`` are the exact lab-frame reconstructions
        (5.1)/(5.3).  Default False reproduces the legacy solver bit-for-bit
        (the alpha=0 backward-compatibility toggle).
    f_r_hz : float, required if displaced
        Bare resonator frequency omega_r [Hz] entering (R.3).
    drive_op_hz : Qobj, required if displaced
        Hz-dimensioned g*(n_hat-n_g) (x) I_r joint operator
        (``charge_drive_operator``); must be Hermitian (asserted).
    """

    def __init__(self, H_cqed, a, kappa_hz, c_ops_extra=None, n_it=3,
                 displaced=False, f_r_hz=None, drive_op_hz=None):
        c_ops = [collapse(kappa_hz, a)]
        if c_ops_extra:
            c_ops = c_ops + list(c_ops_extra)
        self._L0 = np.ascontiguousarray(liouvillian(H_cqed, c_ops).full())
        # Unit-eps drive Liouvillian: Lc(eps_hz) = eps_hz * Lc_unit with
        # Op_t = 2*eps*(a+a_dag); the 2*pi lives in hz_to_qobj (C1).
        Op_unit = hz_to_qobj(2.0) * (a + a.dag())
        self._Lc_unit = np.ascontiguousarray((0.5 * liouvillian(Op_unit)).full())
        from scipy.sparse import csr_matrix
        self._Lc_unit_sp = csr_matrix(self._Lc_unit)
        a_arr = np.asarray(a.full())
        self._a_arr = a_arr
        self._num_arr = a_arr.conj().T @ a_arr
        self._nq, self._nf = (int(d) for d in a.dims[0])
        D = a_arr.shape[0]
        self._D = D
        self._M = D * D
        self._n_it = int(n_it)
        self._displaced = bool(displaced)
        self._kappa_hz = float(kappa_hz)
        self._f_r_hz = None if f_r_hz is None else float(f_r_hz)
        # vec convention: QuTiP column-stacking; vec index k = i + j*D for
        # element (i,j).  Transpose permutation P: k -> (k // D) + (k % D)*D.
        k = np.arange(self._M)
        self._perm = (k // D) + (k % D) * D
        self._diag_idx = np.arange(D) * (D + 1)
        self._trace_row = 0
        # Trace-bordered LU of L0 (unique steady state requires Gamma_r > 0,
        # same uniqueness note as driven_steadystate).
        A0 = self._L0.copy()
        A0[self._trace_row, :] = 0.0
        A0[self._trace_row, self._diag_idx] = 1.0
        self._A0_lu = _sla.lu_factor(A0)
        b0 = np.zeros(self._M, dtype=complex)
        b0[self._trace_row] = 1.0
        self._rho0_nodrive_vec = _sla.lu_solve(self._A0_lu, b0)
        # Verify the GKSL Hermiticity-preservation symmetry L0 = P conj(L0) P
        # on a random vector (enables reusing the +n LUs for the -n blocks).
        rng = np.random.default_rng(12345)
        v = rng.standard_normal(self._M) + 1j * rng.standard_normal(self._M)
        y_direct = self._L0 @ v
        y_sym = np.conj((self._L0 @ np.conj(v[self._perm])))[self._perm]
        sym_err = np.linalg.norm(y_direct - y_sym) / np.linalg.norm(y_direct)
        self._herm_symmetry_ok = bool(sym_err < 1e-10)
        self.herm_symmetry_err = float(sym_err)

        # --- Displaced-frame caching (DERV-02; omega_d-independent) --------
        # L_+- = A_+-(omega_d) * LG with LG = liouvillian(2*pi*g(n_hat-n_g));
        # only the scalar weights change per solve, so LG is cached once.
        self._LG = None
        self._LG_sp = None
        self.lg_symmetry_err = None
        if self._displaced:
            if self._f_r_hz is None or drive_op_hz is None:
                raise ValueError(
                    "displaced=True requires f_r_hz and drive_op_hz "
                    "(charge_drive_operator).")
            _check_drive_op_hermitian(drive_op_hz)
            G_rad = hz_to_qobj(drive_op_hz)    # single 2*pi site (C1)
            self._LG = np.ascontiguousarray(liouvillian(G_rad).full())
            self._LG_sp = csr_matrix(self._LG)
            # Matrix-level App. A identity P conj(LG) P = LG (Hermitian G):
            # with L_- = conj(A_+)*LG this gives L_- = P conj(L_+) P, the
            # exact condition the negative-chain reuse needs.  ASSERTED, not
            # assumed (g=0 zero operator passes trivially).
            yg_direct = self._LG @ v
            yg_sym = np.conj((self._LG @ np.conj(v[self._perm])))[self._perm]
            scale_g = np.linalg.norm(yg_direct)
            lg_err = (np.linalg.norm(yg_direct - yg_sym) / scale_g
                      if scale_g > 0.0 else 0.0)
            self.lg_symmetry_err = float(lg_err)
            if lg_err >= 1e-10:
                raise RuntimeError(
                    "Displaced-mode sideband generator violates "
                    "P conj(LG) P = LG (rel err "
                    f"{lg_err:.2e}); the fast-path negative-chain reuse "
                    "(DERV-02 App. A) is invalid here -- use the production "
                    "driven_steadystate (independent S chain) instead.")

    @property
    def nodrive_n(self) -> float:
        """No-drive <a_dag a> from the cached L0 steady state (diagnostic;
        the formal precondition stays ``nodrive_precondition``)."""
        rho = self._rho0_nodrive_vec.reshape(self._D, self._D, order="F")
        return float(np.real(np.trace(self._num_arr @ rho)))

    def solve(self, omega_d_hz: float, eps_hz: float, n_it: int = None,
              resid_tol: float = 1e-9):
        """Solve the harmonic system at probe omega_d; return sideband/DC data.

        EXACT dense block elimination -- the same downward matrix-continued-
        fraction the production ``driven_steadystate`` runs (T chain
        T_n = -(L0 + i n w + Lc T_{n+1})^{-1} Lc), in dense LAPACK with the
        cached L0 and with the negative-n (S) chain obtained from the GKSL
        Hermiticity-preservation symmetry  S-correction = P conj(Lc T_1) P
        (rho_{-n} = rho_n^dag).  No iterative step; the only differences from
        the production path are dense-vs-CSR linear algebra and the symmetry
        reuse, both cross-validated against ``driven_steadystate`` in the
        Plan-02 check script.

        Returns dict with ``a_wd`` (= <a>(omega_d) = Tr(a rho_+1), feed into
        ``s21``), ``n_dc`` (= Tr(a_dag a rho_0)), ``a_dc`` (~0 self-check),
        and diagnostics ``resid_rel`` (true block-system residual, scaled by
        ||L0 rho_0||) and ``herm_err`` (Hermiticity of rho_0).
        """
        if not self._herm_symmetry_ok:
            raise RuntimeError(
                "GKSL Hermiticity-preservation symmetry check failed at "
                f"construction (rel err {self.herm_symmetry_err:.2e}); the "
                "symmetry-reused negative chain would be invalid.")
        N = self._n_it if n_it is None else int(n_it)
        M = self._M
        w = hz_to_qobj(float(omega_d_hz))      # 2*pi*omega_d [rad/s] (C1 site)
        # Sideband generators for this omega_d (slot mapping per the DERV-02
        # section note: L_plus multiplies e^{-i omega_d t}, L_minus its
        # conjugate harmonic; legacy equal-weight case has L_plus = L_minus).
        if self._displaced:
            alpha_p, alpha_m, A_p = _displaced_sideband_weights(
                self._f_r_hz, float(omega_d_hz), float(eps_hz),
                self._kappa_hz)
            Lp_dense = A_p * self._LG                  # L_plus, dense (RHS)
            Lp_sp = A_p * self._LG_sp                  # L_plus, sparse
            Lm_sp = np.conj(A_p) * self._LG_sp         # L_minus = conj(A_p)*LG
        else:
            alpha_p = 0.0 + 0.0j
            alpha_m = 0.0 + 0.0j
            Lp_dense = float(eps_hz) * self._Lc_unit       # dense (multi-RHS)
            Lp_sp = float(eps_hz) * self._Lc_unit_sp       # sparse (products)
            Lm_sp = Lp_sp                                  # equal weights
        L0 = self._L0
        idx = np.arange(M)
        perm, diag_idx, tr = self._perm, self._diag_idx, self._trace_row

        # Downward T-chain elimination, n = N .. 1 (identical recursion to
        # driven_steadystate; T after the loop is T_1):
        # T_n = -(L0 + i n w + L_minus T_{n+1})^{-1} L_plus.
        lus = []     # [(n, lu_factor(B_n))], n = N .. 1
        T = None
        for n in range(N, 0, -1):
            B = L0.copy()
            B[idx, idx] += 1j * n * w
            if T is not None:
                B += Lm_sp @ T
            lu = _sla.lu_factor(B)
            T = -_sla.lu_solve(lu, Lp_dense)
            lus.append((n, lu))

        # n = 0 effective operator: L0 + L_minus T_1 + (L_plus S_1 by the
        # App. A symmetry: C_minus = P conj(C_plus) P on Hermitian rho_0).
        C_plus = np.asarray(Lm_sp @ T)
        C_minus = np.conj(C_plus[np.ix_(perm, perm)])
        A0 = L0 + C_plus + C_minus
        A0[tr, :] = 0.0
        A0[tr, diag_idx] = 1.0          # trace-bordered row: Tr(rho_0) = 1
        b0 = np.zeros(M, dtype=complex)
        b0[tr] = 1.0
        v0 = _sla.lu_solve(_sla.lu_factor(A0), b0)
        del A0, C_plus, C_minus

        # Chain vectors: rho_n = -B_n^{-1} L_plus rho_{n-1} (n = 1..N);
        # rho_{-n} = rho_n^dag (App. A, valid for V_- = V_+^dag -- asserted
        # at construction in displaced mode).
        vs = {0: v0}
        vprev = v0
        for n, lu in reversed(lus):     # (1, lu_1), (2, lu_2), ..., (N, lu_N)
            vn = -_sla.lu_solve(lu, np.asarray(Lp_sp @ vprev))
            vs[n] = vn
            vprev = vn
        for n in range(1, N + 1):
            vs[-n] = np.conj(vs[n][perm])
        del lus

        # True residual of the truncated harmonic block system (diagnostic);
        # the trace-replaced row of the n=0 block is excluded (its equation
        # was exchanged for Tr(rho_0) = 1, as in steadystate 'direct').
        scale = np.linalg.norm(L0 @ v0)
        rmax = 0.0
        for n in range(-N, N + 1):
            r = L0 @ vs[n] + (1j * n * w) * vs[n]
            if n - 1 >= -N:
                r += Lp_sp @ vs[n - 1]
            if n + 1 <= N:
                r += Lm_sp @ vs[n + 1]
            if n == 0:
                r[tr] = 0.0
            rmax = max(rmax, float(np.linalg.norm(r)))
        resid = rmax / max(scale, 1e-300)
        if not np.isfinite(resid) or resid > resid_tol:
            raise RuntimeError(
                f"FloquetSidebandSolver block elimination residual too large "
                f"at omega_d = {omega_d_hz:.3f} Hz (resid={resid:.2e}); fall "
                "back to driven_steadystate at this point.")

        D = self._D
        rho_0 = v0.reshape(D, D, order="F")
        rho_p1 = vs[1].reshape(D, D, order="F")
        rho_m1 = vs[-1].reshape(D, D, order="F")
        herm_err = float(np.linalg.norm(rho_0 - rho_0.conj().T)
                         / max(np.linalg.norm(rho_0), 1e-300))

        a_res = complex(np.trace(self._a_arr @ rho_p1))
        n_res = float(np.real(np.trace(self._num_arr @ rho_0)))

        # Leakage observables (DERV-02 5.3): p_top over rho_0 and its max
        # over ALL stored harmonics |n| <= N (truncation-error monitor).
        p_top = float(np.real(_top_fock_trace(rho_0, self._nq, self._nf)))
        p_top_max = max(
            abs(_top_fock_trace(vs[n].reshape(D, D, order="F"),
                                self._nq, self._nf))
            for n in range(-N, N + 1))

        if self._displaced:
            # Exact lab-frame reconstruction (DERV-02 (5.1), (5.3)).
            a_wd = alpha_p + a_res
            a_res_m1 = complex(np.trace(self._a_arr @ rho_m1))
            n_dc = (n_res
                    + 2.0 * float(np.real(np.conj(alpha_p) * a_res))
                    + 2.0 * float(np.real(np.conj(alpha_m) * a_res_m1))
                    + abs(alpha_p) ** 2 + abs(alpha_m) ** 2)
        else:
            a_wd = a_res
            n_dc = n_res

        return {
            "a_wd": a_wd,
            "n_dc": n_dc,
            "a_dc": complex(np.trace(self._a_arr @ rho_0)),
            "resid_rel": float(resid),
            "herm_err": herm_err,
            "a_res": a_res,
            "n_res": n_res,
            "p_top": p_top,
            "p_top_max": float(p_top_max),
            "alpha_p": alpha_p,
            "alpha_m": alpha_m,
        }

    def sweep(self, omega_d_grid_hz, eps_hz: float, kappa_c_hz: float,
              n_it: int = None):
        """S21(omega_d) on a grid from the rho_+1 sideband (never DC <a>).

        Same output contract as ``driven_s21_sweep`` (omega_d_hz, a_expect,
        s21, n_expect), fed by the cached-L0 dense solver.
        """
        grid = np.asarray(omega_d_grid_hz, dtype=float)
        a_vals = np.empty(grid.size, dtype=complex)
        n_vals = np.empty(grid.size, dtype=float)
        for j, wd in enumerate(grid):
            r = self.solve(wd, eps_hz, n_it=n_it)
            a_vals[j] = r["a_wd"]
            n_vals[j] = r["n_dc"]
        s21_vals = np.array([s21(av, eps_hz, kappa_c_hz) for av in a_vals])
        return {"omega_d_hz": grid, "a_expect": a_vals, "s21": s21_vals,
                "n_expect": n_vals}


def vald01_peak_scan(solver, f_r_hz, kappa_hz, kappa_c_hz, eps_hz,
                     chi_ref_hz, window_factor=1.6, coarse_points=17,
                     mid_halfwidth_kappa=0.4, mid_step_kappa=0.05,
                     fine_halfwidth_kappa=0.06, fine_step_kappa=0.01,
                     n_it=None):
    """Coarse-locate then densify: staged omega_d scan to f_peak (~0.01*kappa).

    Stage A (sign scan): a coarse grid SYMMETRIC about the bare f_r spanning
    +/- ``window_factor*|chi_ref|`` -- wide enough that a sign-flipped pole
    (the Phase-1 -|chi| failure mode) is FOUND, not assumed away.  |chi_ref|
    (from chi_numeric) sets only the window WIDTH, never the side.
    Stage B: +/- ``mid_halfwidth_kappa``*kappa about the coarse argmax at
    ``mid_step_kappa``*kappa.
    Stage C (densified): +/- ``fine_halfwidth_kappa``*kappa at
    ``fine_step_kappa``*kappa (default ~0.01*kappa, the contract grid), with
    edge re-centering; ``extract_peak_hz`` parabolic sub-grid refinement.

    Returns dict with ``f_peak_hz``, ``n_dc_peak``, the per-stage grids and
    S21 curves, and the +/-|chi_ref| response asymmetry diagnostic.
    """
    W = abs(window_factor * chi_ref_hz)
    coarse_grid = np.linspace(f_r_hz - W, f_r_hz + W, int(coarse_points))
    coarse = solver.sweep(coarse_grid, eps_hz, kappa_c_hz, n_it=n_it)
    resp_c = np.abs(coarse["s21"] - 1.0)
    # Sign-asymmetry diagnostic at +/-|chi_ref| (interpolated).
    resp_plus = float(np.interp(f_r_hz + abs(chi_ref_hz), coarse_grid, resp_c))
    resp_minus = float(np.interp(f_r_hz - abs(chi_ref_hz), coarse_grid, resp_c))
    center_b = float(coarse_grid[int(np.argmax(resp_c))])

    mid_grid = center_b + np.arange(
        -mid_halfwidth_kappa, mid_halfwidth_kappa + 0.5 * mid_step_kappa,
        mid_step_kappa) * kappa_hz
    mid = solver.sweep(mid_grid, eps_hz, kappa_c_hz, n_it=n_it)
    center_c = float(mid_grid[int(np.argmax(np.abs(mid["s21"] - 1.0)))])

    fine = None
    fine_grid = None
    for _ in range(4):  # edge re-centering guard
        fine_grid = center_c + np.arange(
            -fine_halfwidth_kappa, fine_halfwidth_kappa + 0.5 * fine_step_kappa,
            fine_step_kappa) * kappa_hz
        fine = solver.sweep(fine_grid, eps_hz, kappa_c_hz, n_it=n_it)
        k = int(np.argmax(np.abs(fine["s21"] - 1.0)))
        if 0 < k < fine_grid.size - 1:
            break
        center_c = float(fine_grid[k])
    f_peak = extract_peak_hz(fine_grid, fine["s21"])
    k = int(np.argmax(np.abs(fine["s21"] - 1.0)))
    return {
        "f_peak_hz": float(f_peak),
        "n_dc_peak": float(fine["n_expect"][k]),
        "fine_step_hz": float(fine_step_kappa * kappa_hz),
        "asym_ratio": resp_plus / max(resp_minus, 1e-300),
        "coarse": coarse, "mid": mid, "fine": fine,
    }


def vald01_gate_point(qpd, offset_charge, parity, *, coupling_g_hz,
                      readout_freq_hz, kappa_hz, kappa_c_hz, eps_hz,
                      gamma_r_hz, gamma_phi_hz, n_qubit=8, n_fock=6,
                      n_it=3, tol_floor=1e-3, charge_cutoff=30,
                      scan_kwargs=None, return_solver=False):
    """One (n_g, parity) gate point: precondition + densified peak + both legs' inputs.

    Order of operations (the factoring is the point):
      1. ``nodrive_precondition`` (HARD): if it fails, the peak comparison at
         this point is INVALID and is skipped (no peak is extracted from a
         corrupted vacuum).
      2. ``chi_numeric[0]`` from compute_dispersive_matrix(method='numeric',
         rwa=False) on the IDENTICAL build_cqed_hamiltonian parameters
         (same n_g, parity, g, f_r, n_qubit, n_photon = n_fock-1,
         charge_cutoff) -- the leg-(a) anchor.
      3. ``chi_obs[0]``, crossing, weight_frac from
         compute_observed_chi_spectrum on the same parameters -- the leg-(b)
         PR #19 benchmark.
      4. Densified staged peak scan (``vald01_peak_scan``) on the lab-frame
         Floquet sideband S21.

    Residuals (leg a: f_peak - omega_r - chi_numeric[0]; leg b:
    chi_numeric[0] - chi_obs[0]) are returned; tolerance comparison happens in
    ``vald01_run_gate`` once the MEASURED tol is known.
    """
    H, _ = wrap_cqed_hamiltonian(
        qpd, offset_charge, coupling_g_hz, readout_freq_hz, parity=parity,
        n_qubit=n_qubit, n_fock=n_fock, charge_cutoff=charge_cutoff)
    a = cavity_operator(n_qubit, n_fock)
    ops = dressed_jump_operators(
        qpd, offset_charge, coupling_g_hz, readout_freq_hz, parity=parity,
        n_qubit=n_qubit, n_fock=n_fock, charge_cutoff=charge_cutoff)
    c_extra = [collapse(gamma_r_hz, ops["J_minus_dressed"]),
               collapse(gamma_phi_hz / 2.0, ops["Sz_dressed"])]

    # (1) HARD precondition: no-drive vacuum floor (Plan-01 helper).
    nd = nodrive_precondition(H, a, kappa_hz, c_ops_extra=c_extra,
                              tol_floor=tol_floor, assert_floor=False)

    # (2) Leg-(a) anchor on the IDENTICAL Hamiltonian parameters.
    _, chi_ip = qpd.compute_dispersive_matrix(
        offset_charge, coupling_g_hz, readout_freq_hz, num_levels=2,
        parity=parity, charge_cutoff=charge_cutoff, method="numeric",
        n_qubit=n_qubit, n_photon=n_fock - 1, rwa=False)
    chi_numeric0 = float(chi_ip[0])

    # (3) Leg-(b) PR #19 benchmark (same parameters).
    chi_obs, crossing, weight_frac = qpd.compute_observed_chi_spectrum(
        offset_charge, coupling_g_hz, readout_freq_hz, parity=parity,
        n_qubit=n_qubit, n_photon=n_fock - 1, charge_cutoff=charge_cutoff)
    chi_obs0 = float(chi_obs[0])
    crossing0 = bool(crossing[0])
    weight_frac0 = float(weight_frac[0])

    row = {
        "n_g": float(offset_charge), "parity": str(parity),
        "n_qubit": int(n_qubit), "n_fock": int(n_fock), "n_it": int(n_it),
        "nodrive_n": float(nd["n_nd"]), "precondition_passed": bool(nd["passed"]),
        "tol_floor": float(tol_floor),
        "chi_numeric0_hz": chi_numeric0, "chi_obs0_hz": chi_obs0,
        "crossing0": crossing0, "weight_frac0": weight_frac0,
        "f_peak_hz": None, "leg_a_resid_hz": None,
        "leg_b_resid_hz": chi_numeric0 - chi_obs0,
        "sign_ok": None, "asym_ratio": None, "n_dc_peak": None,
        "scan": None,
    }
    solver = None
    if nd["passed"]:
        solver = FloquetSidebandSolver(H, a, kappa_hz, c_ops_extra=c_extra,
                                       n_it=n_it)
        scan = vald01_peak_scan(solver, readout_freq_hz, kappa_hz, kappa_c_hz,
                                eps_hz, chi_numeric0, **(scan_kwargs or {}))
        row["f_peak_hz"] = scan["f_peak_hz"]
        row["leg_a_resid_hz"] = scan["f_peak_hz"] - readout_freq_hz - chi_numeric0
        row["sign_ok"] = bool(np.sign(scan["f_peak_hz"] - readout_freq_hz)
                              == np.sign(chi_numeric0))
        row["asym_ratio"] = scan["asym_ratio"]
        row["n_dc_peak"] = scan["n_dc_peak"]
        row["scan"] = scan
    if return_solver:
        return row, solver
    return row


def _refined_peak(qpd, offset_charge, parity, center_hz, *, coupling_g_hz,
                  readout_freq_hz, kappa_hz, kappa_c_hz, eps_hz, gamma_r_hz,
                  gamma_phi_hz, n_qubit, n_fock, n_it, charge_cutoff=30,
                  fine_step_kappa=0.01):
    """Mid+fine refinement around a known peak at (possibly bumped) settings.

    Used by the measured-tolerance constructor: rebuilds H/dissipators at the
    bumped truncation and re-extracts f_peak around ``center_hz`` (no coarse
    sign scan -- the sign is established by the unbumped staged scan).
    """
    H, _ = wrap_cqed_hamiltonian(
        qpd, offset_charge, coupling_g_hz, readout_freq_hz, parity=parity,
        n_qubit=n_qubit, n_fock=n_fock, charge_cutoff=charge_cutoff)
    a = cavity_operator(n_qubit, n_fock)
    ops = dressed_jump_operators(
        qpd, offset_charge, coupling_g_hz, readout_freq_hz, parity=parity,
        n_qubit=n_qubit, n_fock=n_fock, charge_cutoff=charge_cutoff)
    c_extra = [collapse(gamma_r_hz, ops["J_minus_dressed"]),
               collapse(gamma_phi_hz / 2.0, ops["Sz_dressed"])]
    solver = FloquetSidebandSolver(H, a, kappa_hz, c_ops_extra=c_extra,
                                   n_it=n_it)
    mid_grid = center_hz + np.arange(-0.25, 0.25 + 0.025, 0.05) * kappa_hz
    mid = solver.sweep(mid_grid, eps_hz, kappa_c_hz, n_it=n_it)
    center_c = float(mid_grid[int(np.argmax(np.abs(mid["s21"] - 1.0)))])
    fine = None
    fine_grid = None
    for _ in range(4):
        fine_grid = center_c + np.arange(
            -0.06, 0.06 + 0.5 * fine_step_kappa, fine_step_kappa) * kappa_hz
        fine = solver.sweep(fine_grid, eps_hz, kappa_c_hz, n_it=n_it)
        k = int(np.argmax(np.abs(fine["s21"] - 1.0)))
        if 0 < k < fine_grid.size - 1:
            break
        center_c = float(fine_grid[k])
    return float(extract_peak_hz(fine_grid, fine["s21"]))


def vald01_measured_tolerance(qpd, offset_charge, parity, base_scan, solver,
                              *, coupling_g_hz, readout_freq_hz, kappa_hz,
                              kappa_c_hz, eps_hz, gamma_r_hz, gamma_phi_hz,
                              n_qubit=8, n_fock=6, n_it=3, charge_cutoff=30):
    """MEASURE tol = max(grid_resolution, truncation_drift, stark_term) [Hz].

    All three contributions are measured at a representative (n_g, parity)
    point, NOT pre-asserted (the old 0.05*kappa was finer than the old grid
    -> unmeasurable, the gate's second defect; it is not carried forward):

      * grid_resolution: the densified fine grid step (~0.01*kappa) AND the
        realized sub-grid floor |f_peak(step) - f_peak(step/2)| under a
        halved fine step with the SAME solver -- the larger of the parabolic
        floor and that difference is reported, the step itself is quoted as
        the conservative grid term.
      * truncation_drift: peak shift under n_fock->n_fock+1,
        n_qubit->n_qubit+1, and n_it 3->5 (largest of the three).
      * stark_term: <n>_peak * |dchi/dn| with dchi/dn = chi_0(1) - chi_0(0)
        from compute_stark_spectrum (the AC-Stark/cross-Kerr allowance).
    """
    f0 = float(base_scan["f_peak_hz"])
    step_hz = float(base_scan["fine_step_hz"])

    # --- grid term: halve the fine step around the same peak --------------
    fine_grid = f0 + np.arange(-0.03, 0.03 + 0.0025, 0.005) * kappa_hz
    fine = solver.sweep(fine_grid, eps_hz, kappa_c_hz, n_it=n_it)
    f_half = float(extract_peak_hz(fine_grid, fine["s21"]))
    grid_subgrid_hz = abs(f_half - f0)
    grid_term_hz = step_hz  # conservative: quote the densified step itself

    # --- truncation drifts -------------------------------------------------
    common = dict(coupling_g_hz=coupling_g_hz, readout_freq_hz=readout_freq_hz,
                  kappa_hz=kappa_hz, kappa_c_hz=kappa_c_hz, eps_hz=eps_hz,
                  gamma_r_hz=gamma_r_hz, gamma_phi_hz=gamma_phi_hz,
                  charge_cutoff=charge_cutoff)
    f_fock = _refined_peak(qpd, offset_charge, parity, f0,
                           n_qubit=n_qubit, n_fock=n_fock + 1, n_it=n_it,
                           **common)
    f_nq = _refined_peak(qpd, offset_charge, parity, f0,
                         n_qubit=n_qubit + 1, n_fock=n_fock, n_it=n_it,
                         **common)
    # n_it bump reuses the cached representative solver (same L0).
    fine5 = solver.sweep(f0 + np.arange(-0.06, 0.06 + 0.005, 0.01) * kappa_hz,
                         eps_hz, kappa_c_hz, n_it=5)
    f_nit5 = float(extract_peak_hz(fine5["omega_d_hz"], fine5["s21"]))
    drift_fock_hz = f_fock - f0
    drift_nq_hz = f_nq - f0
    drift_nit_hz = f_nit5 - f0
    truncation_drift_hz = max(abs(drift_fock_hz), abs(drift_nq_hz),
                              abs(drift_nit_hz))

    # --- AC-Stark allowance -------------------------------------------------
    _, chi_n_hz, _ = qpd.compute_stark_spectrum(
        offset_charge, coupling_g_hz, readout_freq_hz, num_levels=2,
        parity=parity, n_qubit=n_qubit, n_photon=n_fock - 1,
        charge_cutoff=charge_cutoff)
    dchi_dn_hz = float(chi_n_hz[0, 1] - chi_n_hz[0, 0])
    n_peak = float(base_scan["n_dc_peak"])
    stark_term_hz = abs(n_peak * dchi_dn_hz)

    tol_hz = max(grid_term_hz, truncation_drift_hz, stark_term_hz)
    return {
        "grid_step_hz": step_hz,
        "grid_subgrid_floor_hz": grid_subgrid_hz,
        "grid_term_hz": grid_term_hz,
        "drift_fock_hz": drift_fock_hz,
        "drift_nqubit_hz": drift_nq_hz,
        "drift_nit_hz": drift_nit_hz,
        "truncation_drift_hz": truncation_drift_hz,
        "dchi_dn_hz_per_photon": dchi_dn_hz,
        "n_dc_peak": n_peak,
        "stark_term_hz": stark_term_hz,
        "tol_hz": tol_hz,
        "tol_a_hz": tol_hz,
        "tol_b_hz": tol_hz,
        "measurable": bool(tol_hz >= step_hz),
    }


def vald01_run_gate(qpd, *, coupling_g_hz, readout_freq_hz, kappa_hz,
                    kappa_c_hz, eps_hz, gamma_r_hz, gamma_phi_hz,
                    n_g_points=(0.0, 0.25, 0.5), parities=("even", "odd"),
                    n_qubit=8, n_fock=6, n_it=3, tol_floor=1e-3,
                    tol_point=(0.0, "even"), charge_cutoff=30, verbose=True):
    """Run the tightened, factored VALD-01 gate across the (n_g, parity) set.

    Per point: hard no-drive precondition -> densified staged peak scan ->
    leg-(a) residual vs chi_numeric on the identical H -> leg-(b) residual vs
    PR #19 chi_obs (crossing/weight_frac asserted).  The tolerance is MEASURED
    at ``tol_point`` and applied to all points (tol_a = tol_b = measured max).

    Returns dict: ``rows`` (per-point records incl. scan curves),
    ``tolerance`` (measured breakdown), and the verdict booleans.  A FAIL is
    a valid result and is reported, never tuned away (fp-magnitude).
    """
    import time as _time

    common = dict(coupling_g_hz=coupling_g_hz, readout_freq_hz=readout_freq_hz,
                  kappa_hz=kappa_hz, kappa_c_hz=kappa_c_hz, eps_hz=eps_hz,
                  gamma_r_hz=gamma_r_hz, gamma_phi_hz=gamma_phi_hz,
                  n_qubit=n_qubit, n_fock=n_fock, n_it=n_it,
                  tol_floor=tol_floor, charge_cutoff=charge_cutoff)
    rows = []
    rep_solver = None
    rep_scan = None
    for n_g in n_g_points:
        for parity in parities:
            t0 = _time.perf_counter()
            row, solver = vald01_gate_point(qpd, n_g, parity,
                                            return_solver=True, **common)
            row["wall_s"] = _time.perf_counter() - t0
            rows.append(row)
            if (n_g, parity) == tuple(tol_point):
                rep_solver = solver
                rep_scan = row["scan"]
            if verbose:
                fp = (f"{row['f_peak_hz'] - readout_freq_hz:+.1f}"
                      if row["f_peak_hz"] is not None else "INVALID")
                print(f"  point n_g={n_g} {parity}: nodrive <n>="
                      f"{row['nodrive_n']:.3e} ({'ok' if row['precondition_passed'] else 'FAIL'})"
                      f"  f_peak-f_r={fp} Hz  chi_num={row['chi_numeric0_hz']:+.1f}"
                      f"  chi_obs={row['chi_obs0_hz']:+.1f}  [{row['wall_s']:.0f} s]")

    if rep_solver is None or rep_scan is None:
        raise RuntimeError(
            f"tolerance point {tol_point} has no valid scan (precondition "
            "failed?); cannot measure tol -- gate INVALID, do not assert a pass.")

    tol = vald01_measured_tolerance(
        qpd, tol_point[0], tol_point[1], rep_scan, rep_solver,
        coupling_g_hz=coupling_g_hz, readout_freq_hz=readout_freq_hz,
        kappa_hz=kappa_hz, kappa_c_hz=kappa_c_hz, eps_hz=eps_hz,
        gamma_r_hz=gamma_r_hz, gamma_phi_hz=gamma_phi_hz,
        n_qubit=n_qubit, n_fock=n_fock, n_it=n_it, charge_cutoff=charge_cutoff)

    tol_a = tol["tol_a_hz"]
    tol_b = tol["tol_b_hz"]
    for row in rows:
        valid = row["precondition_passed"]
        row["leg_a_pass"] = bool(valid and row["leg_a_resid_hz"] is not None
                                 and abs(row["leg_a_resid_hz"]) < tol_a
                                 and row["sign_ok"])
        row["leg_b_pass"] = bool((not row["crossing0"])
                                 and row["weight_frac0"] >= 0.9
                                 and abs(row["leg_b_resid_hz"]) < tol_b)
    precondition_all = all(r["precondition_passed"] for r in rows)
    leg_a_all = all(r["leg_a_pass"] for r in rows)
    leg_b_all = all(r["leg_b_pass"] for r in rows)
    sign_all = all(bool(r["sign_ok"]) for r in rows)
    gate_pass = precondition_all and leg_a_all and leg_b_all and sign_all
    return {
        "rows": rows, "tolerance": tol, "tol_a_hz": tol_a, "tol_b_hz": tol_b,
        "precondition_all": precondition_all, "leg_a_all": leg_a_all,
        "leg_b_all": leg_b_all, "sign_all": sign_all, "gate_pass": gate_pass,
        "rep_solver": rep_solver, "rep_scan": rep_scan,
    }


def vald01_eps_extrapolation(solver, f_r_hz, kappa_hz, kappa_c_hz,
                             eps_list_hz, center_hz, chi_obs0_hz,
                             fine_halfwidth_kappa=0.08, fine_step_kappa=0.01,
                             n_it=None):
    """eps->0 extrapolation: f_peak(<n>) -> omega_r + chi_obs[0] (self-evidencing).

    Runs the densified fine scan at each drive amplitude in ``eps_list_hz`` on
    the SAME absolute omega_d grid (so the parabolic sub-grid bias is common
    mode and cancels in differences), records (<n>_peak, f_peak), then fits
    f_peak = intercept + slope*<n>.  The intercept is the linear-response
    pole (must land on omega_r + chi_obs[0] within tol_a); the slope is the
    AC-Stark/cross-Kerr pull (small, negative; ~ dchi/dn).
    """
    grid = center_hz + np.arange(
        -fine_halfwidth_kappa, fine_halfwidth_kappa + 0.5 * fine_step_kappa,
        fine_step_kappa) * kappa_hz
    eps_arr = np.asarray(eps_list_hz, dtype=float)
    f_peaks = np.empty(eps_arr.size)
    n_peaks = np.empty(eps_arr.size)
    curves = []
    for j, eps in enumerate(eps_arr):
        sw = solver.sweep(grid, eps, kappa_c_hz, n_it=n_it)
        f_peaks[j] = extract_peak_hz(grid, sw["s21"])
        k = int(np.argmax(np.abs(sw["s21"] - 1.0)))
        n_peaks[j] = sw["n_expect"][k]
        curves.append(sw)
    slope, intercept = np.polyfit(n_peaks, f_peaks, 1)
    target = f_r_hz + chi_obs0_hz
    return {
        "eps_hz": eps_arr, "n_dc_peak": n_peaks, "f_peak_hz": f_peaks,
        "intercept_hz": float(intercept), "slope_hz_per_photon": float(slope),
        "target_hz": float(target),
        "intercept_minus_target_hz": float(intercept - target),
        "grid_hz": grid, "curves": curves,
    }


# ===========================================================================
# Phase 3 Plan 03-01: RESEARCH-licensed light extensions (ADDITIVE ONLY)
# ===========================================================================
#
# Three drivers/analysis routines AROUND the locked Hamiltonian -- NO new
# physics, NO re-derivation of DERV-01/02, and NOTHING above touched.  Every
# routine consumes the existing full-cosine multilevel charge-basis pieces
# (``wrap_cqed_hamiltonian`` / ``charge_drive_operator`` /
# ``dressed_jump_operators`` / ``cavity_operator`` / ``s21``) verbatim:
#
#   (1) ``mesolve_branch_run``        -- long-time qutip.mesolve of the SAME
#       lab-frame H_lab(t) = H_cQED + 2*eps*cos(omega_d t)(a+a_dag) from a
#       chosen initial state (vacuum or a bright/displaced coherent state).
#       The steady state of a single-photon (linear) cavity drive is UNIQUE
#       (Bartolo et al., PRA 94, 033841 (2016)); hysteresis can therefore only
#       arise from IC-dependent TRANSIENT dynamics, which is exactly what this
#       engine exposes (vacuum-branch vs bright-branch plateaus).
#
#   (2) ``semiclassical_floquet_onset`` -- branch-analysis / quasienergy
#       onset on the EXISTING charge basis (Dumas-Cohen-...-Blais, PRX 14,
#       041023 (2024), arXiv:2402.06615): a semiclassical mean-field Duffing
#       bracket for the bifurcation scale PLUS non-secular Floquet quasienergy
#       crossings (built from the SAME steadystate_fourier harmonic structure,
#       never QuTiP's fmmesolve/FloquetBasis) to flag the multiphoton-resonance
#       onset photon number.  Ionization n_crit != JC n_crit.
#
#   (3) ``cavity_qfunc`` / ``pointer_snr`` -- thin wrappers over qutip.qfunc
#       for pointer-lobe / bimodality readout and the parity pointer SNR.
#
# HARD CONSTRAINTS (the binding Phase-3 locks):
#   * fp-fmmesolve (driven.py:698-708): NO fmmesolve / FloquetBasis /
#     floquet_tensor / steadystate_floquet anywhere -- the secular step is a
#     hidden RWA-like reduction violating the non-RWA lock.  The onset routine
#     uses the in-core non-secular Fourier recursion ONLY.
#   * fp-twolevel: every routine runs on the full-cosine multilevel charge
#     basis (no nq=2 hardcode, no g^2/Delta dispersive substitution, no RWA
#     drive); E_C/anharmonicity sensitivity is preserved (test-no-twolevel).
#   * The 02-03-validated machinery (driven_steadystate, FloquetSidebandSolver,
#     displacement_alphas, charge_drive_operator, dressed_jump_operators, s21)
#     is UNCHANGED -- these routines only call it, never edit it.  The 02-03
#     permanent regression fixture is the gate on this file (Task 2).
# ===========================================================================


def mesolve_branch_run(
    H_cqed,
    a,
    c_ops,
    psi0,
    tlist_s,
    eps_hz: float,
    omega_d_hz: float,
    *,
    options=None,
):
    """Long-time lab-frame mesolve from a chosen IC (vacuum/bright branch engine).

    Integrates the SAME lab-frame master equation the steady-state solver
    encodes,

        H_lab(t) = H_cQED + 2*eps*cos(omega_d t)*(a + a_dag),

    with ``qutip.mesolve`` from the initial state ``psi0`` (a ket or density
    matrix) over ``tlist_s`` and reads off the long-time plateaus
    ``<n>`` = Tr(a_dag a rho(t_end)) and ``<a>`` = Tr(a rho(t_end)).  Running
    this from a VACUUM IC and from a BRIGHT (displaced/coherent) IC and
    comparing the plateaus is the bistability/hysteresis discriminator: the
    single-photon-drive steady state is UNIQUE (Bartolo 2016), so below the
    onset BOTH ICs must reach the SAME plateau (no spurious bistability), while
    IC-dependent plateaus at finite power flag a long-lived metastable branch.

    This is ADDITIVE: it does NOT touch the validated displaced/Floquet
    steady-state machinery -- it is an independent dynamical cross-check on the
    same locked H.  There is NO -omega_d*a_dag a rotating-frame shift
    (fp-rf-shortcut): the coupling stays static+exact in the lab frame and the
    only time dependence is the cosine cavity drive (the factor 2 makes each
    e^{-+i omega_d t} sideband carry the C3-locked eps*(a+a_dag), matching the
    steady-state convention bit-for-bit).

    All frequency inputs are in Hz; the single 2*pi site is ``hz_to_qobj``
    (drive amplitude and time-domain phase alike) and the collapse rates are
    already 2*pi-scaled inside ``c_ops`` (built via ``collapse``).  The cosine
    time argument uses omega_d in rad/s = 2*pi*omega_d_hz (the same single
    conversion); ``tlist_s`` is in SECONDS so that 2*pi*f*t is dimensionless.

    Parameters
    ----------
    H_cqed : Qobj
        Static joint Hamiltonian in rad/s (``wrap_cqed_hamiltonian``), used AS
        IS -- lab frame, non-RWA, full cosine; NO omega_d shift.
    a : Qobj
        Bare cavity annihilation in the joint space (``cavity_operator``).
    c_ops : list of Qobj
        Collapse operators, already 2*pi-scaled (``collapse(kappa, a)`` plus
        the C5 dressed set).  At least one qubit relaxation channel is needed
        for relaxation to a unique long-time state.
    psi0 : Qobj
        Initial state (ket or density matrix) in the joint space -- e.g.
        ``tensor(basis(nq,0), basis(nf,0))`` (vacuum) or a bright coherent
        cavity state ``tensor(basis(nq,0), coherent(nf, alpha))``.
    tlist_s : array_like
        Time grid in SECONDS.  Must extend to several cavity lifetimes
        (>> 1/(2*pi*kappa)) for the plateau to be the steady state.
    eps_hz : float
        Drive amplitude epsilon in Hz (C3 normalization; lab cos drive is
        2*eps).
    omega_d_hz : float
        Probe/drive frequency omega_d in Hz.
    options : dict or qutip.Options, optional
        Passed to ``qutip.mesolve`` (e.g. tighter atol/rtol).

    Returns
    -------
    dict with keys
        ``n_plateau`` : float, long-time <a_dag a> at t = tlist_s[-1];
        ``a_plateau`` : complex, long-time <a> at t = tlist_s[-1];
        ``n_traj`` : ndarray, <a_dag a>(t) over tlist_s;
        ``a_traj`` : ndarray (complex), <a>(t) over tlist_s;
        ``tlist_s`` : ndarray, the time grid echoed back;
        ``rho_final`` : Qobj, the final-time state (for Q-function / pointer
                        analysis).
    """
    from qutip import mesolve

    num = a.dag() * a
    w_d_rad = hz_to_qobj(float(omega_d_hz))     # 2*pi*omega_d [rad/s], C1 site
    # Lab-frame drive coefficient: 2*eps*cos(omega_d t), the 2*pi already in
    # both the amplitude (hz_to_qobj) and the cosine argument (w_d_rad).
    drive_amp_rad = hz_to_qobj(2.0 * float(eps_hz))   # 2*pi*2*eps [rad/s]
    drive_op = a + a.dag()

    # COMPILED string coefficient (QuTiP 5): cos(w*t) with w = omega_d in
    # rad/s passed via args.  A pure-Python callback f(t) is invoked at EVERY
    # integrator substep and dominates the runtime on the joint Liouvillian;
    # the string form compiles to C and is ~100x faster.  The single 2*pi
    # site is preserved (w = 2*pi*omega_d_hz from hz_to_qobj).
    H_t = [H_cqed, [drive_amp_rad * drive_op, "cos(w*t)"]]
    args = {"w": float(w_d_rad)}

    # Lab-frame time-domain integration carries the drive at omega_d, so the
    # integrator must resolve the carrier PERIOD T_d = 2*pi/w_d_rad.  We give
    # zvode a generous nsteps budget and a max_step at the Nyquist bound
    # (half the drive period) so the cosine is never aliased, with QuTiP's
    # DEFAULT tolerances (atol 1e-8 / rtol 1e-6) -- tighter tolerances force
    # pathologically small adaptive steps on the joint Liouvillian for no
    # physics gain (Deviation Rule 2, numerical: these are integrator-accuracy
    # controls, NOT physics changes; the master equation, lab frame, and
    # convention are untouched).  Caller-supplied ``options`` win.  For a
    # realistic GHz carrier over many cavity lifetimes this lab-frame engine
    # is intentionally a CROSS-CHECK on the validated steady-state solver, not
    # the production sweep path; reduced-carrier limit checks keep it cheap.
    drive_period_s = 2.0 * np.pi / float(w_d_rad)
    default_opts = {
        "nsteps": 100000,
        "max_step": 0.5 * drive_period_s,
    }
    if options:
        default_opts.update(dict(options))

    res = mesolve(H_t, psi0, np.asarray(tlist_s, dtype=float),
                  c_ops=list(c_ops), e_ops=[num, a], args=args,
                  options=default_opts)
    n_traj = np.real(np.asarray(res.expect[0], dtype=float))
    a_traj = np.asarray(res.expect[1], dtype=complex)
    rho_final = res.final_state
    # LAB-FRAME PHASE NOTE.  <a>(t) carries the drive carrier e^{-i omega_d t},
    # so the instantaneous complex a_plateau = <a>(t_end) is the steady
    # amplitude times an arbitrary carrier phase e^{-i omega_d t_end}.  The
    # PHASE-INVARIANT branch discriminants are <n> (= Tr(a_dag a rho), exactly
    # the steady-state occupation, no carrier) and |<a>| (the demodulated
    # amplitude |alpha_+| of the rho_+1 sideband the steady-state solver
    # returns).  Bistability/hysteresis is read from n_plateau / abs_a_plateau,
    # never from the raw complex phase.  ``a_demod`` removes the carrier using
    # the final time so it can be compared to the steady-state -2i eps/kappa.
    tl_arr = np.asarray(tlist_s, dtype=float)
    t_end = float(tl_arr[-1])
    a_demod = complex(a_traj[-1] * np.exp(1j * float(w_d_rad) * t_end))
    # TAIL-MEAN plateau (carrier-ripple-robust).  The INSTANTANEOUS <n>(t_end)
    # still carries a small residual ripple at the drive harmonics from the
    # finite carrier resolution; averaging <n>(t) over the last ``tail_frac``
    # of the window gives the steady occupation cleanly (the demodulated <a>
    # ripple averages out in <n> = Tr(a_dag a rho), which has no carrier).
    tail_frac = 0.2
    n_tail = max(2, int(tail_frac * n_traj.size))
    n_plateau_mean = float(np.mean(n_traj[-n_tail:]))
    return {
        "n_plateau": float(n_traj[-1]),
        "n_plateau_mean": n_plateau_mean,
        "a_plateau": complex(a_traj[-1]),
        "abs_a_plateau": float(abs(a_traj[-1])),
        "a_demod": a_demod,
        "n_traj": n_traj,
        "a_traj": a_traj,
        "tlist_s": tl_arr,
        "rho_final": rho_final,
    }


def semiclassical_floquet_onset(
    qpd,
    offset_charge: float,
    coupling_g_hz: float,
    readout_freq_hz: float,
    kappa_hz: float,
    *,
    parity: str = "even",
    n_qubit: int = 8,
    n_fock: int = 13,
    charge_cutoff: int = 30,
    n_levels_floquet: int = 6,
    eps_probe_hz: float = None,
    omega_d_hz: float = None,
    n_it: int = 5,
):
    """Branch-analysis + non-secular Floquet quasienergy onset (Blais 2024 method).

    Runs on the EXISTING full-cosine multilevel charge basis (no two-level /
    Kerr / RWA reduction; fp-twolevel) and combines two onset estimators:

    (A) SEMICLASSICAL DUFFING BRACKET.  From the undriven dressed spectrum
        (``solve_cqed_eigensystem``) extract the qubit anharmonicity and the
        dispersive (cross-Kerr) cavity pull, giving a mean-field Duffing/Kerr
        coefficient ``K`` for the driven cavity mode.  The classical Duffing
        bistability threshold sets a bracketing photon scale
        ``n_bracket ~ kappa / (sqrt(3) * |K|)`` -- an ORDER-OF-MAGNITUDE bracket
        on the bifurcation, NOT the JC n_crit (ionization onset differs; Blais).

    (B) NON-SECULAR FLOQUET QUASIENERGY CROSSINGS.  Build the in-core
        steadystate_fourier harmonic structure at a weak probe (the SAME
        non-secular recursion the validated solver uses -- NEVER QuTiP's
        fmmesolve/FloquetBasis/floquet_tensor, fp-fmmesolve) and locate the
        multiphoton-resonance onset by where dressed transition frequencies
        cross integer multiples of omega_d (the lab-frame quasienergy folding).
        Returns the lowest photon number at which a |0,n> level enters a
        k-photon resonance with a higher charge-basis level.

    The two estimators bracket the finite-power onset from the classical
    (Duffing) and quantum (level-crossing) sides; neither is the JC n_crit and
    both are reported as ESTIMATES (the decisive onset is a Phase 3-02/03
    deliverable, not this infra plan).  E_C sensitivity is structural (the
    anharmonicity and the level spacings both move with E_C); the
    ``test-no-twolevel`` check exercises exactly that.

    Parameters
    ----------
    qpd : qpd.theory.transmon.QPD
        Transmon instance (E_J, E_C in Hz) -- the full-cosine CPB.
    offset_charge, parity
        Gate point (C6).
    coupling_g_hz, readout_freq_hz, kappa_hz : float
        Device coupling g, bare resonator f_r, total cavity kappa [Hz].
    n_qubit, n_fock, charge_cutoff : int
        Truncation (full multilevel charge basis; no reduction).
    n_levels_floquet : int
        Number of low-lying transmon levels scanned for k-photon crossings.
    eps_probe_hz, omega_d_hz : float, optional
        If both given, also evaluate the in-core non-secular Floquet harmonic
        response at this probe (cross-check the recursion runs on the same H);
        defaults skip the explicit solve and use the spectral crossing analysis
        only.
    n_it : int
        Harmonic truncation for any in-core Floquet evaluation.

    Returns
    -------
    dict with keys
        ``anharmonicity_hz`` : float, qubit 12-01 anharmonicity (E_C scale);
        ``cross_kerr_hz`` : float, dispersive cavity pull dchi/dn proxy;
        ``kerr_coeff_hz`` : float, mean-field Duffing K;
        ``n_bracket_duffing`` : float, Duffing bistability photon bracket;
        ``floquet_onset_photons`` : float, lowest k-photon-resonance photon
                                    number from the level-crossing scan
                                    (inf if no crossing within the scanned
                                    levels);
        ``crossing_table`` : list of (level, k_photons, n_photons) records;
        ``basis_dim`` : int, joint Hilbert dimension actually used (proves the
                        multilevel charge basis was consumed);
        ``used_fmmesolve`` : bool, always False (fp-fmmesolve guard echo).
    """
    # --- Undriven dressed spectrum on the FULL charge basis (no reduction) --
    evals_hz, evecs, labels, _ = qpd.solve_cqed_eigensystem(
        offset_charge, coupling_g_hz, readout_freq_hz, parity=parity,
        n_qubit=n_qubit, n_photon=n_fock - 1, rwa=False,
        charge_cutoff=charge_cutoff)
    evals_hz = np.asarray(evals_hz, dtype=float)
    labels = list(labels)
    basis_dim = int(evals_hz.size)

    # Index the |qubit=i, photon=0> ladder from the dressed labels.
    def _level0(i):
        for k, lab in enumerate(labels):
            if lab[0] == i and lab[1] == 0:
                return k
        return None

    i0, i1, i2 = _level0(0), _level0(1), _level0(2)
    # Qubit 0->1 and 1->2 transition frequencies (undriven, dressed).
    f01 = float(evals_hz[i1] - evals_hz[i0]) if (i1 is not None) else float("nan")
    f12 = (float(evals_hz[i2] - evals_hz[i1])
           if (i2 is not None and i1 is not None) else float("nan"))
    anharmonicity_hz = f12 - f01    # ~ -E_C for a transmon (E_C-sensitive)

    # Dispersive cavity pull (cross-Kerr) from the observed-chi spectrum slope
    # dchi/dn (the AC-Stark per-photon shift; the validated chi machinery).
    _, chi_n_hz, _ = qpd.compute_stark_spectrum(
        offset_charge, coupling_g_hz, readout_freq_hz, num_levels=2,
        parity=parity, n_qubit=n_qubit, n_photon=n_fock - 1,
        charge_cutoff=charge_cutoff)
    cross_kerr_hz = float(chi_n_hz[0, 1] - chi_n_hz[0, 0])   # dchi/dn

    # Mean-field Duffing/Kerr coefficient: the cavity self-Kerr induced by the
    # qubit nonlinearity is ~ the per-photon cross-Kerr pull (the cavity
    # frequency moves by ~cross_kerr per photon).  Duffing bistability needs
    # the nonlinear shift over the population to exceed the linewidth:
    # n_bracket ~ kappa / (sqrt(3) |K|) (classical Duffing threshold).
    kerr_coeff_hz = cross_kerr_hz
    if abs(kerr_coeff_hz) > 0.0:
        n_bracket = float(kappa_hz / (np.sqrt(3.0) * abs(kerr_coeff_hz)))
    else:
        n_bracket = float("inf")

    # --- (B) k-photon level-crossing scan (lab-frame quasienergy folding) ---
    # A |0,n> photon-dressed level enters a k-photon resonance with a higher
    # charge-basis level j when E_j - E_0 ~ k * (omega_r + n*cross_kerr) for
    # integer k.  Scan the low-lying transmon levels and report the lowest
    # photon number n at which such a crossing first occurs.
    crossing_table = []
    onset_photons = float("inf")
    omega_r = float(readout_freq_hz)
    n_scan = np.arange(0, 4 * max(int(np.ceil(abs(n_bracket)
                                              if np.isfinite(n_bracket) else 8)),
                                  8) + 1)
    for j in range(1, min(n_levels_floquet, basis_dim)):
        ij = _level0(j)
        if ij is None:
            continue
        dE = float(evals_hz[ij] - evals_hz[i0])  # j-level excitation energy
        # Effective drive photon energy at occupation n (AC-Stark-shifted).
        for n_ph in n_scan:
            drive_energy = omega_r + n_ph * cross_kerr_hz
            if drive_energy <= 0:
                continue
            k_real = dE / drive_energy
            k_round = int(round(k_real))
            if k_round >= 1 and abs(k_real - k_round) < 0.02:
                crossing_table.append((int(j), int(k_round), int(n_ph)))
                if n_ph < onset_photons:
                    onset_photons = float(n_ph)
                break

    # --- Optional in-core non-secular Floquet evaluation (recursion echo) ---
    # Proves the SAME non-secular steadystate_fourier path runs on this H; it
    # is NOT QuTiP's fmmesolve (forbidden).  Skipped unless a probe is given.
    floquet_probe = None
    if eps_probe_hz is not None and omega_d_hz is not None:
        H, comps = wrap_cqed_hamiltonian(
            qpd, offset_charge, coupling_g_hz, readout_freq_hz, parity=parity,
            n_qubit=n_qubit, n_fock=n_fock, charge_cutoff=charge_cutoff)
        a = cavity_operator(n_qubit, n_fock)
        ops = dressed_jump_operators(
            qpd, offset_charge, coupling_g_hz, readout_freq_hz, parity=parity,
            n_qubit=n_qubit, n_fock=n_fock, charge_cutoff=charge_cutoff)
        gamma_r = kappa_hz / 100.0
        c_extra = [collapse(gamma_r, ops["J_minus_dressed"]),
                   collapse(gamma_r / 2.0, ops["Sz_dressed"])]
        res = driven_steadystate(
            H, a, float(omega_d_hz), float(eps_probe_hz), kappa_hz,
            c_ops_extra=c_extra, n_it=n_it)
        floquet_probe = {"a_wd": res["a_wd"], "n_dc": res["n_dc"],
                         "p_top_max": res["p_top_max"]}

    return {
        "anharmonicity_hz": float(anharmonicity_hz),
        "f01_hz": float(f01),
        "f12_hz": float(f12),
        "cross_kerr_hz": float(cross_kerr_hz),
        "kerr_coeff_hz": float(kerr_coeff_hz),
        "n_bracket_duffing": float(n_bracket),
        "floquet_onset_photons": float(onset_photons),
        "crossing_table": crossing_table,
        "floquet_probe": floquet_probe,
        "basis_dim": basis_dim,
        "used_fmmesolve": False,
    }


# ===========================================================================
# DERV-03 (plan 03.1-01): ADDITIVE non-secular Sambe/Shirley Floquet solver.
#
# This block is an ADDITIVE extension only. It introduces NO change to the
# existing semiclassical_floquet_onset arm above, nor to any Phase-2/3
# validated path (driven_steadystate, FloquetSidebandSolver, displacement_
# alphas, charge_drive_operator, dressed_jump_operators, s21). It adds a
# self-contained, hand-rolled Sambe/extended-space (Shirley) quasienergy
# engine for the SEMICLASSICAL classical-cavity Floquet arm of Phase 3.1.
#
# Convention source of truth (pinned, confirmed 2026-06-17):
#   (pinned Floquet convention chain)
#   CONVENTIONS-floquet.md
#
#   * Drive on the TRANSMON: eps_d[Hz] = g[Hz]*sqrt(nbar_r), acting on the
#     (n_hat - n_g_eff) charge operator in the truncated transmon eigenbasis
#     (the n_shift_mat carried by build_cqed_hamiltonian return_components --
#     NOT the bare n_hat from _qubit_block_in_eigenbasis, and NOT the cavity
#     AC-Stark pull eps_cavity=(kappa/2)*sqrt(<n>); fp-eps-confusion). The
#     semiclassical arm has NO quantized cavity -- the field is replaced by
#     its classical mean amplitude alpha(t)=sqrt(nbar)*e^{-i omega_d t}.
#   * Drive time-dependence (CONVENTIONS sec 3, Blais/Dumas 2024 form):
#       H(t) = H_t - i*eps_d*sin(omega_d t)*(n_hat - n_g_eff)
#     Fourier-expanding -i*sin(omega_d t) = -(1/2)(e^{+i omega_d t}
#     - e^{-i omega_d t}). The two k=+/-1 Sambe off-diagonal blocks are
#     therefore IMAGINARY: H_F[m, m+1] = -i*(eps_d/2)*(n_hat - n_g_eff),
#     H_F[m, m-1] = (H_F[m+1, m])^dag, making H_F manifestly Hermitian for a
#     Hermitian (real-symmetric) charge matrix. CONVENTIONS sec 7 writes the
#     blocks in the schematic +/-(eps_d/2)N form (as for a 2V*cos drive); for
#     the literal -i*sin drive and a real-symmetric N that schematic real form
#     is NOT Hermitian, so the implementation uses the manifestly-Hermitian
#     -i*(eps_d/2)N blocks. This is a phase-origin (sin-vs-cos) gauge choice
#     ONLY: the quasienergy spectrum is byte-identical to the real-block
#     gauge to ~1e-11 relative (verified), per CONVENTIONS sec 3.
#   * Non-secular by construction (fp-fmmesolve): a single dense
#     scipy.linalg.eigh of the assembled Sambe matrix. NO fmmesolve /
#     FloquetBasis / floquet_tensor / steadystate_floquet anywhere; no
#     time-averaging / secular step. used_fmmesolve=False is echoed.
#   * Full-cosine multilevel charge basis (fp-twolevel): the transmon block
#     comes from the repo CPB primitive via solve_eigensystem; no two-level /
#     Kerr / fixed-E_J/E_C reduction. basis_dim is echoed.
# ===========================================================================
def sambe_floquet_quasienergies(
    qpd,
    nbar_r: float,
    coupling_g_hz: float,
    omega_d_hz: float,
    *,
    offset_charge: float = 0.0,
    parity: str = "even",
    n_qubit: int = 8,
    charge_cutoff: int = 30,
    n_harmonics: int = 10,
    eps_d_hz: float = None,
):
    """Non-secular Sambe/Shirley Floquet quasienergies of the driven transmon.

    Hand-rolled extended-space (Sambe/Shirley) solver for the semiclassical
    classically-driven full-cosine QPD transmon (DERV-03, plan 03.1-01). The
    drive ``eps_d = g*sqrt(nbar_r)`` acts DIRECTLY on the transmon charge
    operator ``(n_hat - n_g_eff)`` so the transmon-ionization Stark shift is
    present NON-secularly (the knob the 03-01 quasienergy heuristic missed).
    The convention chain is pinned in CONVENTIONS-floquet.md (confirmed
    2026-06-17); this routine is the source-of-truth implementation of its
    section (7) Sambe matrix.

    The semiclassical arm has no quantized cavity: the resonator field is
    replaced by its classical mean amplitude ``alpha(t)=sqrt(nbar)*
    e^{-i omega_d t}`` so the only Hilbert space is the truncated transmon
    eigenbasis. The Sambe matrix is built on
    ``{transmon} (x) {Fourier modes m in [-n_harmonics, n_harmonics]}``:

        H_F[m, m]    =  H_t + m*omega_d*I
        H_F[m, m+1]  =  -i*(eps_d/2)*(n_hat - n_g_eff)
        H_F[m, m-1]  =  (H_F[m+1, m])^dag      (Hermitian conjugate block)

    diagonalized once with dense ``scipy.linalg.eigh``. The blocks are
    imaginary because the drive is ``-i*eps_d*sin(omega_d t)`` (CONVENTIONS
    sec 3); the spectrum is gauge-identical to the real-block (cos phase
    origin) form.

    Parameters
    ----------
    qpd : qpd.theory.transmon.QPD
        Transmon instance (E_J, E_C in Hz) -- the full-cosine CPB. Use the
        CALC-01/Phase-3 device point E_J=6.95 GHz, E_C=0.695 GHz; NOT the
        CLAUDE.md E_J=8.335 GHz example.
    nbar_r : float
        Cavity mean photon number (the drive-strength variable). Sets
        ``eps_d = g*sqrt(nbar_r)`` unless ``eps_d_hz`` is given explicitly.
        This is the SAME nbar used in the 314 mean-field bracket (no hidden
        sqrt(2) or <n>-vs-|alpha|^2 factor; CONVENTIONS sec 5).
    coupling_g_hz : float
        Bare device coupling g [Hz].
    omega_d_hz : float
        Drive frequency [Hz] (the readout drive f_r at the device point).
    offset_charge : float
        Dimensionless offset charge n_g (period 1 = 2e). The parity shift
        (+0.5 for 'odd') is applied inside the CPB primitive (C6); the
        ``n_g_eff`` subtraction carried on the charge operator is
        ``offset_charge + (0.5 if parity=='odd' else 0.0)``.
    parity : {'even', 'odd'}
        CPB parity. Primary gate point is ``parity='even'``, ``offset_charge=0``
        (n_g_eff=0), apples-to-apples with the 314 bracket and CALC-01.
    n_qubit : int
        Number of transmon eigenstates kept (>= ~8 so above-well / ionized
        crossing partners are represented; convergence knob for plan 03.1-02).
    charge_cutoff : int
        CPB charge-basis cutoff passed to the primitive (full-cosine basis).
    n_harmonics : int
        Sambe Fourier-mode half-width ``N_h``: modes ``m in [-N_h, N_h]``
        (convergence knob for plan 03.1-02). Sambe matrix size is
        ``n_qubit*(2*N_h+1)`` -- a trivial dense ``eigh``.
    eps_d_hz : float, optional
        Override the transmon drive amplitude directly [Hz]. If given,
        ``nbar_r`` is ignored for the drive (still echoed). Used by the
        Stark-off limiting case (``eps_d_hz=0.0``).

    Returns
    -------
    dict with keys
        ``quasienergies_hz`` : ndarray, shape (n_qubit*(2*N_h+1),)
            Sorted Sambe eigenvalues [Hz] (NOT folded -- the caller folds
            mod omega_d as needed).
        ``modes`` : ndarray, complex, shape (n_qubit*(2*N_h+1),
            2*N_h+1, n_qubit)
            For each eigenvector j, ``modes[j]`` is the eigenvector reshaped
            to (harmonic m, transmon level) so branch-overlap following is
            possible. Row index ``m_index = m + N_h``.
        ``m0_weight`` : ndarray, shape (n_qubit*(2*N_h+1),)
            Total probability weight in the m=0 Fourier block for each mode
            (used to identify the physical/computational branch).
        ``eps_d_hz`` : float, the transmon drive amplitude used [Hz].
        ``nbar_r`` : float, the photon number used.
        ``omega_d_hz`` : float.
        ``qubit_freqs_hz`` : ndarray, the bare transmon eigenfrequencies [Hz].
        ``n_shift_mat`` : ndarray, the (n_hat - n_g_eff) charge matrix used.
        ``n_harmonics`` : int, ``N_h``.
        ``basis_dim`` : int, ``n_qubit*(2*N_h+1)`` (the full Sambe dimension;
            proves the multilevel charge basis was consumed, no reduction).
        ``n_g_eff`` : float, the offset actually subtracted.
        ``used_fmmesolve`` : bool, always False (fp-fmmesolve guard echo).
    """
    if nbar_r < 0.0:
        raise ValueError(f"nbar_r must be >= 0, got {nbar_r}")
    n_g_eff = float(offset_charge) + (0.5 if parity == "odd" else 0.0)

    # --- Transmon block + (n_hat - n_g_eff) on the FULL-COSINE charge basis ---
    # Reuse the repo CPB primitive (no two-level / Kerr reduction; fp-twolevel).
    # _qubit_block_in_eigenbasis applies the parity shift internally and returns
    # the BARE n_hat in the eigenbasis; we apply the n_g_eff subtraction here to
    # form the SAME (n_hat - n_g_eff) operator build_cqed_hamiltonian uses
    # (CONVENTIONS sec 4). At n_g_eff=0 (even) the subtraction is a no-op.
    qubit_freqs_hz, n_bare_mat = qpd._qubit_block_in_eigenbasis(
        offset_charge, parity=parity, n_qubit=n_qubit, charge_cutoff=charge_cutoff
    )
    qubit_freqs_hz = np.asarray(qubit_freqs_hz, dtype=float)
    n_shift_mat = np.asarray(n_bare_mat, dtype=float) - n_g_eff * np.eye(n_qubit)

    # The off-diagonal Stark coupling lives in the off-diagonal matrix elements
    # of n_shift_mat, which the diagonal n_g_eff subtraction never touches; the
    # subtraction only adds a state-independent diagonal piece (CONVENTIONS sec 4).
    if not np.allclose(n_shift_mat, n_shift_mat.T, atol=1e-10):
        raise ValueError(
            "charge operator (n_hat - n_g_eff) is not symmetric to 1e-10; "
            "the -i*sin Sambe blocks require a Hermitian charge matrix"
        )

    eps_d = (float(eps_d_hz) if eps_d_hz is not None
             else float(coupling_g_hz) * np.sqrt(float(nbar_r)))

    # --- Assemble the Sambe/Shirley matrix (Hz units; single 2*pi site, C1,
    # is irrelevant here because every term -- H_t, m*omega_d, eps_d -- carries
    # the SAME frequency units, so the eigenvalues come out directly in Hz and
    # multiplying the whole matrix by 2*pi would only rescale them uniformly). -
    N_h = int(n_harmonics)
    dim = int(n_qubit)
    n_modes = 2 * N_h + 1
    big = n_modes * dim
    H_t = np.diag(qubit_freqs_hz).astype(complex)
    # H_F[m, m+1] = -i*(eps_d/2)*(n_hat - n_g_eff)  [coeff of e^{-i omega_d t}
    # in the -i*sin expansion gives the +1 super-diagonal block]
    V_up = -1j * (eps_d / 2.0) * n_shift_mat.astype(complex)
    H_F = np.zeros((big, big), dtype=complex)
    eye_dim = np.eye(dim)
    for a, m in enumerate(range(-N_h, N_h + 1)):
        sl = slice(a * dim, (a + 1) * dim)
        H_F[sl, sl] = H_t + m * float(omega_d_hz) * eye_dim
        if a + 1 < n_modes:
            su = slice((a + 1) * dim, (a + 2) * dim)
            H_F[sl, su] = V_up            # block (m, m+1)
            H_F[su, sl] = V_up.conj().T   # block (m+1, m) = Hermitian conjugate

    # Hermiticity guard (the -i*sin blocks must combine Hermitian; eigvals real).
    herm_resid = np.max(np.abs(H_F - H_F.conj().T))
    if herm_resid > 1e-6:
        raise ValueError(
            f"Sambe matrix not Hermitian (resid={herm_resid:.3e} Hz); "
            "check the drive-block convention"
        )

    # --- Single dense diagonalization (non-secular by construction). ---
    # NO fmmesolve / FloquetBasis / floquet_tensor / steadystate_floquet; no
    # time-averaging / secular step. This IS the non-secular spectrum.
    evals, evecs = _sla.eigh(H_F)  # dense, real eigenvalues (Hermitian)
    order = np.argsort(evals)
    evals = evals[order]
    evecs = evecs[:, order]

    # Reshape each eigenvector to (harmonic m, transmon level) for branch follow.
    modes = np.empty((big, n_modes, dim), dtype=complex)
    for j in range(big):
        modes[j] = evecs[:, j].reshape(n_modes, dim)
    m0_index = N_h
    m0_weight = np.sum(np.abs(modes[:, m0_index, :]) ** 2, axis=1)

    return {
        "quasienergies_hz": evals,
        "modes": modes,
        "m0_weight": m0_weight,
        "eps_d_hz": float(eps_d),
        "nbar_r": float(nbar_r),
        "omega_d_hz": float(omega_d_hz),
        "qubit_freqs_hz": qubit_freqs_hz,
        "n_shift_mat": n_shift_mat,
        "n_harmonics": N_h,
        "basis_dim": int(big),
        "n_g_eff": float(n_g_eff),
        "used_fmmesolve": False,
    }


def sambe_floquet_branch_scan(
    qpd,
    nbar_grid,
    coupling_g_hz: float,
    omega_d_hz: float,
    *,
    comp_level: int = 0,
    offset_charge: float = 0.0,
    parity: str = "even",
    n_qubit: int = 8,
    charge_cutoff: int = 30,
    n_harmonics: int = 10,
    stark_off: bool = False,
    p_comp_crossing: float = 0.5,
):
    """Follow the computational Floquet branch across a photon-number grid.

    Thin branch-follower on top of :func:`sambe_floquet_quasienergies`
    (DERV-03 branch composition). At the FIRST (smallest) nbar in
    ``nbar_grid`` the computational branch is identified as the Floquet mode
    whose m=0 Fourier block is dominantly the bare computational transmon
    state ``|comp_level>`` (max overlap with the m=0 unit vector
    e_{comp_level}). The branch is then tracked across the grid by maximum
    state overlap of the m=0 block (Shillito max-overlap following), and
    ``P_comp(nbar) = |<comp | Phi_branch>|^2`` is recorded -- the projection of
    the tracked branch's m=0 block onto the bare computational state.

    Parameters
    ----------
    qpd, coupling_g_hz, omega_d_hz, offset_charge, parity, n_qubit,
    charge_cutoff, n_harmonics
        Forwarded to :func:`sambe_floquet_quasienergies`.
    nbar_grid : array-like
        Photon numbers to scan (ascending; the first is the branch seed).
    comp_level : int
        Bare transmon level defining the computational branch (0 = ground).
    stark_off : bool
        If True, force ``eps_d=0`` at every grid point (the Stark-off /
        03-01 cavity-AC-Stark-only limit, test-stark-off). The transmon
        levels then never Stark-shift, so no above-well crossing folds in.
    p_comp_crossing : float
        Threshold on P_comp below which the tracked branch is flagged as
        having undergone an avoided crossing (left the computational manifold).

    Returns
    -------
    dict with keys
        ``nbar_grid`` : ndarray, the scanned photon numbers.
        ``eps_d_hz`` : ndarray, eps_d at each grid point.
        ``p_comp`` : ndarray, P_comp(nbar) of the tracked branch.
        ``branch_quasienergy_hz`` : ndarray, the tracked branch quasienergy
            [Hz] at each grid point (the eigenvalue, NOT folded).
        ``partner_level`` : ndarray of int, the dominant transmon level of the
            tracked branch's m=0 block at each grid point (the crossing partner
            once P_comp drops).
        ``crossing_nbar`` : float, the first nbar where P_comp < p_comp_crossing
            (inf if no crossing within the grid -- the Stark-off signature).
        ``comp_level`` : int.
        ``stark_off`` : bool.
        ``used_fmmesolve`` : bool, always False.
    """
    nbar_grid = np.asarray(nbar_grid, dtype=float)
    if nbar_grid.ndim != 1 or nbar_grid.size == 0:
        raise ValueError("nbar_grid must be a non-empty 1-D array")

    eps_d_list = []
    p_comp_list = []
    branch_qe_list = []
    partner_list = []
    crossing_nbar = float("inf")
    prev_m0 = None  # tracked branch m=0 block (complex vector, len n_qubit)

    for nbar in nbar_grid:
        res = sambe_floquet_quasienergies(
            qpd, float(nbar), coupling_g_hz, omega_d_hz,
            offset_charge=offset_charge, parity=parity, n_qubit=n_qubit,
            charge_cutoff=charge_cutoff, n_harmonics=n_harmonics,
            eps_d_hz=(0.0 if stark_off else None),
        )
        modes = res["modes"]          # (big, n_modes, n_qubit)
        N_h = res["n_harmonics"]
        m0_blocks = modes[:, N_h, :]  # (big, n_qubit) -- each mode's m=0 block

        if prev_m0 is None:
            # Seed: pick the mode whose m=0 block is most |comp_level>-like AND
            # carries dominant total m=0 weight (a physical, near-static branch).
            comp_vec = np.zeros(n_qubit, dtype=complex)
            comp_vec[comp_level] = 1.0
            overlaps = np.abs(m0_blocks @ comp_vec.conj()) ** 2
            # Require the branch to be m=0-dominant (physical Floquet state).
            score = overlaps * res["m0_weight"]
            jbr = int(np.argmax(score))
        else:
            # Track by max overlap of the m=0 block with the previous branch.
            overlaps = np.abs(m0_blocks @ prev_m0.conj()) ** 2
            jbr = int(np.argmax(overlaps))

        m0 = m0_blocks[jbr]
        m0_norm = m0 / (np.linalg.norm(m0) + 1e-300)
        # P_comp = projection of the (m=0-normalized) branch onto |comp_level>.
        p_comp = float(np.abs(m0_norm[comp_level]) ** 2)
        partner = int(np.argmax(np.abs(m0_norm) ** 2))

        eps_d_list.append(res["eps_d_hz"])
        p_comp_list.append(p_comp)
        branch_qe_list.append(float(res["quasienergies_hz"][jbr]))
        partner_list.append(partner)
        if p_comp < p_comp_crossing and not np.isfinite(crossing_nbar):
            crossing_nbar = float(nbar)
        prev_m0 = m0_norm

    return {
        "nbar_grid": nbar_grid,
        "eps_d_hz": np.asarray(eps_d_list, dtype=float),
        "p_comp": np.asarray(p_comp_list, dtype=float),
        "branch_quasienergy_hz": np.asarray(branch_qe_list, dtype=float),
        "partner_level": np.asarray(partner_list, dtype=int),
        "crossing_nbar": crossing_nbar,
        "comp_level": int(comp_level),
        "stark_off": bool(stark_off),
        "used_fmmesolve": False,
    }


def cavity_qfunc(rho, xvec, yvec, *, n_qubit=None, n_fock=None):
    """Thin wrapper over ``qutip.qfunc`` for the cavity Husimi-Q (pointer lobes).

    If ``rho`` is a JOINT qubit+cavity state, the qubit is traced out first
    (``rho.ptrace(1)`` with the cavity as subsystem 1, qubit-major C1 order) so
    the Q-function is the cavity reduced state -- the pointer/bimodality
    readout observable.  A single-lobed Q-function (one maximum) is a coherent
    pointer; two lobes signal a parity/branch pointer-state split.

    Parameters
    ----------
    rho : Qobj
        Cavity density matrix, OR a joint qubit+cavity density matrix (then the
        qubit is traced out).
    xvec, yvec : array_like
        Phase-space grids (Re/Im of alpha) passed to ``qutip.qfunc``.
    n_qubit, n_fock : int, optional
        If given and ``rho`` is joint, used only as a dimensional sanity check.

    Returns
    -------
    dict with keys
        ``Q`` : 2D ndarray, the Husimi-Q distribution on (xvec, yvec);
        ``xvec``, ``yvec`` : the grids echoed back;
        ``n_lobes`` : int, number of local maxima above 25% of the global max
                      (1 = coherent pointer; 2 = bimodal branch split).
    """
    from qutip import qfunc

    rho_c = rho
    # Joint state? qubit-major dims [[nq, nf], ...]; trace out the qubit (sub 0)
    # leaving the cavity (sub 1).
    if isinstance(rho.dims[0], list) and len(rho.dims[0]) == 2:
        rho_c = rho.ptrace(1)

    Q = np.asarray(qfunc(rho_c, np.asarray(xvec, dtype=float),
                         np.asarray(yvec, dtype=float)))

    # Count local maxima above 25% of the peak (cheap bimodality proxy).
    thresh = 0.25 * float(Q.max())
    n_lobes = 0
    ny, nx = Q.shape
    for iy in range(1, ny - 1):
        for ix in range(1, nx - 1):
            v = Q[iy, ix]
            if v < thresh:
                continue
            nb = Q[iy - 1:iy + 2, ix - 1:ix + 2]
            if v >= nb.max() - 1e-15:
                n_lobes += 1
    return {
        "Q": Q,
        "xvec": np.asarray(xvec, dtype=float),
        "yvec": np.asarray(yvec, dtype=float),
        "n_lobes": int(max(n_lobes, 1)),
    }


def pointer_snr(alpha_e: complex, alpha_o: complex, field_variance: float):
    """Parity-pointer separation SNR = |alpha_e - alpha_o| / sqrt(variance).

    The dispersive readout distinguishes the even/odd charge parity by the
    cavity pointer position; the single-shot SNR is the pointer separation in
    units of the pointer width.  ``field_variance`` is the per-quadrature
    pointer variance (vacuum = 1/2 in the |alpha> normalization here); the
    caller supplies the measured/steady-state value.

    Parameters
    ----------
    alpha_e, alpha_o : complex
        Even- and odd-parity pointer centers (e.g. <a> from the two parity
        steady states, or the two Q-function lobe centers).
    field_variance : float
        Pointer variance (per quadrature; > 0).

    Returns
    -------
    dict with keys
        ``snr`` : float, |alpha_e - alpha_o| / sqrt(field_variance);
        ``separation`` : float, |alpha_e - alpha_o|;
        ``field_variance`` : float, echoed back.
    """
    sep = float(abs(complex(alpha_e) - complex(alpha_o)))
    var = float(field_variance)
    if var <= 0.0:
        raise ValueError(f"field_variance must be > 0, got {var}")
    return {
        "snr": float(sep / np.sqrt(var)),
        "separation": sep,
        "field_variance": var,
    }
