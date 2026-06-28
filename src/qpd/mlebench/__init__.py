"""MLE-bench dataset infrastructure for the QPD VNA simulator.

Generates two competition-style datasets of complex S21 (I + iQ) time traces
from :mod:`qpd.simulator`:

* **signal**     -- background telegraph + charge events **+ quasiparticle bursts**
* **background** -- background telegraph + charge events, **no bursts**

Both datasets share the same background condition (electronic white noise,
quiescent quasiparticle tunneling, and charge events) and are generated at
equal length in 30 s chunks. The reconstruction targets are:

1. the timing of every parity flip (quasiparticle tunneling),
2. the time and fractional part of each charge jump, and
3. the onset time and quasiparticle count of each burst (signal only).

See :mod:`qpd.mlebench.config` for the frozen physics parameters,
:mod:`qpd.mlebench.generate` for the generator, and
:mod:`qpd.mlebench.grade` for the submission grader.
"""

from qpd.mlebench.config import (
    DEFAULT_PHYSICS,
    DatasetConfig,
    GenerationBudget,
    PhysicsConfig,
)
from qpd.mlebench.generate import (
    ChunkTruth,
    generate_chunk,
    generate_dataset,
    plan_generation,
)
from qpd.mlebench.grade import grade, grade_report

__all__ = [
    "PhysicsConfig",
    "DEFAULT_PHYSICS",
    "DatasetConfig",
    "GenerationBudget",
    "ChunkTruth",
    "generate_chunk",
    "generate_dataset",
    "plan_generation",
    "grade",
    "grade_report",
]
