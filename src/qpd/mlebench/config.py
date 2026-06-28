"""Frozen physics + dataset configuration for the MLE-bench QPD datasets.

The physics numbers are copied verbatim from ``notebooks/simulation.ipynb``
(QPD device, notch resonator, readout LO, noise, telegraph rates and burst EMG
profile). Changing any value here changes the generated dataset, so the whole
config is captured into each dataset's ``metadata.json`` for provenance.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field


@dataclass(frozen=True)
class PhysicsConfig:
    """All physics/simulator parameters, matching the simulation notebook.

    Times are seconds, frequencies/rates are Hz, energies are Hz.
    """

    # --- QPD device (notebook s1) ---
    e_j_hz: float = 8.335e9
    e_c_hz: float = 0.695e9
    coupling_g_hz: float = 150e6

    # --- Notch resonator fit (notebook s2) ---
    fr_dressed_hz: float = 6489353036.703264  # measured DRESSED even-parity line
    q_i: float = 39669.674
    q_c_abs: float = 81560.897
    q_l: float = 26688.765  # informational; recomputed from q_i, q_c_abs
    phi: float = 0.03
    a: float = 0.5
    alpha: float = 0.4
    tau_env: float = 50e-9
    # Bare-cavity back-solve (fixed point on chi): iterations and the +offset
    # used in the notebook's loop ``F_BARE = FR - chi + 15e6``.
    bare_solve_iters: int = 6
    bare_solve_offset_hz: float = 15e6

    # --- Readout LO (notebook s4): parked 1 kHz above the dressed line ---
    f_drive_offset_hz: float = 1e3

    # --- Acquisition ---
    sample_rate_hz: float = 1e5

    # --- Background quasiparticle telegraph (notebook s5) ---
    gamma_even_to_odd_hz: float = 100.0
    gamma_odd_to_even_hz: float = 100.0

    # --- Electronic white noise (per-quadrature std, electronic coords) ---
    noise_sigma: float = 2e-4

    # --- Slow offset-charge drift: sawtooth (notebook s3) ---
    sawtooth_n_g_min: float = -0.5
    sawtooth_n_g_max: float = 5.0
    sawtooth_slope_per_s: float = 50.0

    # --- Charge-jump events (shared background condition) ---
    # Poisson arrival rate of charge jumps; per-jump size is the FRACTIONAL
    # part only -- the integer part of an offset-charge jump is physically
    # unobservable (the CPB Hamiltonian is periodic in n_g with period 1), so
    # the reconstruction target is the fractional jump drawn here.
    charge_jump_rate_hz: float = 0.4
    charge_jump_frac_low: float = -0.5  # uniform draw on (low, high]
    charge_jump_frac_high: float = 0.5

    # --- Quasiparticle bursts (signal dataset only), EMG profile (notebook s5)
    burst_tau: float = 3.7e-3
    burst_mu: float = 1.2e-3
    burst_sigma: float = 0.4e-3
    burst_expected_n_qp: float = 15.0
    # Burst-onset Poisson rate. The user spec is a mean of ~30 onsets per 30 s
    # chunk -> 1.0 Hz. (Background dataset overrides this to 0.)
    burst_onset_rate_hz: float = 1.0

    # --- chi(n_g) interpolation grid (passed to VNASimulator) ---
    n_g_grid_points: int = 401
    charge_cutoff: int = 30

    def f_drive_hz(self) -> float:
        return self.fr_dressed_hz + self.f_drive_offset_hz

    def solve_bare_fr_hz(self, qpd) -> float:
        """Back-solve the bare cavity so the dressed even line lands on fr.

        Mirrors the fixed-point loop in notebook s2.
        """
        f_bare = self.fr_dressed_hz
        for _ in range(self.bare_solve_iters):
            _, chi = qpd.compute_dispersive_matrix(
                0.0,
                self.coupling_g_hz,
                f_bare,
                num_levels=2,
                parity="even",
            )
            f_bare = self.fr_dressed_hz - chi[0] + self.bare_solve_offset_hz
        return float(f_bare)

    def to_dict(self) -> dict:
        return asdict(self)


DEFAULT_PHYSICS = PhysicsConfig()


@dataclass
class DatasetConfig:
    """Per-dataset switches.

    Two datasets are produced from the same :class:`PhysicsConfig`:
    ``signal`` (bursts on) and ``background`` (bursts off). Everything else --
    noise, telegraph, charge jumps -- is identical.
    """

    name: str  # "signal" | "background"
    bursts_enabled: bool
    chunk_seconds: float = 30.0
    train_fraction: float = 0.6  # fraction of chunks whose labels are public

    @property
    def competition_id(self) -> str:
        return f"qpd-{self.name}"


def default_dataset_configs(
    chunk_seconds: float = 30.0, train_fraction: float = 0.6
) -> list[DatasetConfig]:
    return [
        DatasetConfig(
            name="signal",
            bursts_enabled=True,
            chunk_seconds=chunk_seconds,
            train_fraction=train_fraction,
        ),
        DatasetConfig(
            name="background",
            bursts_enabled=False,
            chunk_seconds=chunk_seconds,
            train_fraction=train_fraction,
        ),
    ]


@dataclass
class GenerationBudget:
    """Stop condition for a full (two-dataset) generation run.

    The job stops at whichever limit is reached first: total bytes written or
    total wall-clock compute. The two datasets are kept at EQUAL chunk count,
    so each gets half of the total budget when planning.
    """

    max_total_bytes: float = 2e9  # 2 GB
    max_total_compute_seconds: float = 3600.0  # 1 hour
    # Hard cap on chunks per dataset (safety net; None = budget-driven only).
    max_chunks_per_dataset: int | None = None

    def to_dict(self) -> dict:
        return {
            "max_total_bytes": self.max_total_bytes,
            "max_total_compute_seconds": self.max_total_compute_seconds,
            "max_chunks_per_dataset": self.max_chunks_per_dataset,
        }
