"""Blind reconstruction of charge-parity dynamics from readout traces.

There are two entry points, and which one applies is set by how the offset
charge is driven during the measurement:

:func:`reconstruct_parity_flips_ramped`
    ``n_g`` is **swept** (the sawtooth of ``notebooks/simulation.ipynb``). The
    two branch means move with the ramp, sweep through parity-blind crossings,
    and are relabelled at every sawtooth reset. The routine learns the fold
    period, the reset comb and any offset-charge jumps from the trace.

:func:`reconstruct_parity_flips_static`
    ``n_g`` is **held constant**. The readout collapses to two stationary blobs
    in the I/Q plane with telegraph switching between them, so the entire ramp
    apparatus is unnecessary: fit the two blobs, decode, done. Much simpler and
    much cheaper -- but the contrast is fixed by the chosen bias point instead
    of sweeping, so a bias near the parity-blind charge is unrecoverable.

Both are blind: no device, resonator, or noise parameter is supplied, and both
share the two-state HMM in :mod:`.hmm`. Parity is taken to be a pure telegraph
process -- it holds until the next quasiparticle tunnels, and each tunnel
toggles it -- with equal dwell times for the two states.

:mod:`.analysis` scores a reconstruction against the truth, including the
per-burst efficiency and timing bias as a function of how many tunnels land in
the window.
"""

from qpd.reconstruction.analysis import (
    FlipScore,
    burst_report,
    score_flips,
    window_report,
)
from qpd.reconstruction.benchmark import (
    BenchmarkReport,
    BenchmarkTrial,
    BurstSizePoint,
    TraceFidelity,
    as_complex_trace,
    benchmark_reconstruction,
    benchmark_vs_noise,
    characterize_trace,
    sweep_burst_size,
    sweep_rate,
)
from qpd.reconstruction.benchmark_plots import (
    plot_burst_efficiency,
    plot_burst_multiplicity,
    plot_efficiency_vs_rate,
)
from qpd.reconstruction.bursts import (
    BurstMatch,
    DetectedBurst,
    detect_bursts,
    match_bursts,
)
from qpd.reconstruction.emission import (
    EmissionModel,
    estimate_direction,
    estimate_fold_period,
    learn_emission_model,
)
from qpd.reconstruction.events import extract_flips, flip_confidence
from qpd.reconstruction.plotting import plot_iq_plane, plot_trace_with_flips
from qpd.reconstruction.hmm import (
    HMMResult,
    decode,
    decode_with_rate,
    decoded_path,
    forward_backward,
    viterbi,
)
from qpd.reconstruction.ramp import ResetComb, find_reset_comb
from qpd.reconstruction.reconstruct import (
    ReconstructionResult,
    reconstruct_parity_flips_ramped,
)
from qpd.reconstruction.segment import misfit_statistic, segment_and_realign
from qpd.reconstruction.static_bias import (
    StaticBlobModel,
    StaticReconstructionResult,
    fit_two_blobs,
    reconstruct_parity_flips_static,
)

__all__ = [
    # swept n_g
    "reconstruct_parity_flips_ramped",
    "ReconstructionResult",
    "EmissionModel",
    "learn_emission_model",
    "estimate_direction",
    "estimate_fold_period",
    "segment_and_realign",
    "misfit_statistic",
    "ResetComb",
    "find_reset_comb",
    # constant n_g
    "reconstruct_parity_flips_static",
    "StaticReconstructionResult",
    "StaticBlobModel",
    "fit_two_blobs",
    # shared
    "extract_flips",
    "flip_confidence",
    "HMMResult",
    "decode",
    "decode_with_rate",
    "decoded_path",
    "forward_backward",
    "viterbi",
    # plotting
    "plot_trace_with_flips",
    "plot_iq_plane",
    "FlipScore",
    "score_flips",
    "window_report",
    "burst_report",
    # benchmarking a measured trace
    "characterize_trace",
    "TraceFidelity",
    "benchmark_reconstruction",
    "benchmark_vs_noise",
    "BenchmarkReport",
    "BenchmarkTrial",
    "as_complex_trace",
    # burst finding and burst-level diagnostics
    "detect_bursts",
    "match_bursts",
    "DetectedBurst",
    "BurstMatch",
    "sweep_rate",
    "sweep_burst_size",
    "BurstSizePoint",
    "plot_efficiency_vs_rate",
    "plot_burst_efficiency",
    "plot_burst_multiplicity",
]
