"""Blind reconstruction of charge-parity dynamics from readout traces.

The entry point is :func:`reconstruct_parity_flips`, which recovers the timing
of every parity flip from a complex I/Q trace without being told any device,
resonator, or noise parameter -- the discrimination axis, the noise level, the
two-branch splitting profile, the offset-charge ramp and its reset schedule are
all learned from the trace.

Parity is taken to be a pure telegraph process: it holds until the next
quasiparticle tunnels, and each tunnel toggles it. Both dwell times are drawn
from the same rate, so the model is symmetric.

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
from qpd.reconstruction.emission import (
    EmissionModel,
    estimate_direction,
    estimate_fold_period,
    learn_emission_model,
)
from qpd.reconstruction.events import extract_flips, flip_confidence
from qpd.reconstruction.hmm import HMMResult, decode, forward_backward, viterbi
from qpd.reconstruction.ramp import ResetComb, find_reset_comb
from qpd.reconstruction.reconstruct import (
    ReconstructionResult,
    reconstruct_parity_flips,
)
from qpd.reconstruction.segment import misfit_statistic, segment_and_realign

__all__ = [
    "reconstruct_parity_flips",
    "ReconstructionResult",
    "EmissionModel",
    "learn_emission_model",
    "estimate_direction",
    "estimate_fold_period",
    "segment_and_realign",
    "misfit_statistic",
    "ResetComb",
    "find_reset_comb",
    "extract_flips",
    "flip_confidence",
    "HMMResult",
    "decode",
    "forward_backward",
    "viterbi",
    "FlipScore",
    "score_flips",
    "window_report",
    "burst_report",
]
