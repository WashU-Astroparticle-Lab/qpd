# qpd
Theoretical calculations for Quantum Parity Detectors.

## Install

```bash
pip install -e .
```

## What is here

- **`qpd.theory`** — the QPD transmon model. Energy levels, dispersive shifts,
  quantum-capacitance fitting, materials database.
  See [`docs/theory.md`](docs/theory.md) and [`docs/materials.md`](docs/materials.md).
- **`qpd.simulator`** — forward model: generates realistic I/Q readout traces
  with telegraph parity switching, quasiparticle bursts, offset-charge drift and
  jumps, and a Probst notch resonator response.
- **`qpd.reconstruction`** — the inverse problem: recovers the timing of every
  quasiparticle tunnelling event from an I/Q trace, and benchmarks its own
  efficiency and accuracy on a measured trace by surrogate replay.
  See [`docs/reconstruction.md`](docs/reconstruction.md).
- **`qpd.mlebench`** — dataset generator and grader built on the above.

```python
from qpd import QPD

qpd = QPD(e_j_hz=8.335e9, e_c_hz=0.695e9)
qpd.plot_all()
```

## Reconstruction

Given a readout trace, recover when the parity flipped. Blind by default: the
discrimination axis, the noise level, the two-branch splitting and the
offset-charge ramp are all learned from the trace, with no device or resonator
parameter supplied.

```python
from qpd.reconstruction import (reconstruct_parity_flips_static,
                                reconstruct_parity_flips_ramped,
                                plot_trace_with_flips, score_flips)

# n_g held constant -> two stationary blobs in the I/Q plane
rec = reconstruct_parity_flips_static(iq, sample_rate=1e5)

# n_g swept (sawtooth) -> moving branches, blind crossings, ramp resets
rec = reconstruct_parity_flips_ramped(iq, sample_rate=1e5)

rec.flip_times        # reconstructed tunnelling times [s]
rec.degenerate        # True => output should not be trusted; see the notes
rec.decoded_fidelity  # expected per-sample assignment fidelity
```

Always check `rec.degenerate` and `rec.contrast` before using the output: a
model that has latched onto noise fails *quietly*, and the fidelity estimate
stays high while it does. [`docs/reconstruction.md`](docs/reconstruction.md)
§13b explains what to check and why.

### How well does it work on *your* data?

Measured data has no truth, so `score_flips` cannot run on it.
`benchmark_reconstruction` fits your trace, replays that fitted fidelity into
surrogate traces that *do* have truth, and reconstructs those blind with the
same settings — so the efficiency and accuracy it reports are those of the data
you actually took.

```python
from qpd.reconstruction import benchmark_reconstruction

report = benchmark_reconstruction(iq, sample_rate=1e5, n_trials=16)
print(report.summary())

report.efficiency          # (mean, sd) recall
report.purity              # (mean, sd) precision
report.timing_rms_s        # (mean, sd)
report.corrected_rate_hz   # measured flip count, de-biased by the two above
report.warnings            # read these first
```

Read `report.warnings` before quoting any of it: on a degenerate trace the
fitted model is spurious, so the surrogates replay a *fiction* that
reconstructs well.
[`docs/reconstruction.md`](docs/reconstruction.md) §12c has the method and the
evidence that it predicts.

### Rate response and burst response

Three diagnostics, all driven by the fidelity of the same measured trace.

```python
from qpd.reconstruction import (characterize_trace, sweep_rate, sweep_burst_size,
                                plot_efficiency_vs_rate, plot_burst_efficiency,
                                plot_burst_multiplicity)

fid = characterize_trace(iq, 1e5)

# 1. flip efficiency vs Poisson background rate
plot_efficiency_vs_rate(sweep_rate(fid, [1, 10, 100, 1e3, 1e4]))

# 2, 3. burst-LEVEL detection efficiency and multiplicity bias vs burst size,
#       at the background rate measured from the data
points = sweep_burst_size(fid, [2, 3, 5, 8, 12, 20, 50], background_rate_hz=bg)
plot_burst_efficiency(points)
plot_burst_multiplicity(points)
```

`detect_bursts` clusters a reconstructed flip train against the Poisson
background with a trials-corrected scan statistic (measured false-burst rate
< 0.01 per 5 s trace). Note that recovered burst multiplicity **saturates** —
crowded tunnels cancel in pairs — so a large burst's quoted count is a lower
bound. [`docs/reconstruction.md`](docs/reconstruction.md) §12d has the figures.

## Verification

```bash
python checks/check_parity_reconstruction.py     # reconstruction regression gate
python checks/check_reconstruction_benchmark.py  # surrogate-replay benchmark gate
python checks/check_readout_window.py            # dressed-S21 + lock-in window
python checks/study_parity_reconstruction.py     # rate / noise / burst / device sweeps
python checks/study_reconstruction_diagnostics.py  # the three figures above
```

## Notebooks

`notebooks/simulation.ipynb` runs the simulator end to end and, in §7,
demonstrates the reconstruction against the simulator's own ground truth.
