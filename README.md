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

## LED-pulsed waveforms

`LEDPulseModel` generates synchronized, independent chunks: the first pulse is
at 50 ms, followed by 50 µs ON and 300 ms OFF. Each pulse draws a Poisson
tunneling count and reuses the existing Gaussian-plus-exponential time profile,
with excitation spread over ON and delays constrained to be nonnegative.
Offset charge is redrawn uniformly on `[-0.5, 0.5)` at each rising edge and at
independent random jumps averaging one per minute. Quiescent tunneling remains
active at 20 Hz throughout ON and OFF.

```python
from qpd import QPD, ResonatorConfig, VNASimulator, LEDPulseModel, led_count_sweep

device = QPD(e_j_hz=6.95e9, e_c_hz=0.695e9)
device.coupling_g_hz = 40e6
sim = VNASimulator(
    qpd=device, resonator=ResonatorConfig(f_r=6.936e9, q_i=1e5, q_c_abs=1.2e5),
    f_drive=6.936437e9, sample_rate=1e5,
    led_pulsing=LEDPulseModel(expected_n_qp=10),
)
result = sim.simulate(duration=5.0, seed=0)
count_settings = led_count_sweep()  # zero + 31 log-spaced means from 1 to 1000
```

`expected_n_qp` counts **tunneling events per pulse**, not physical particles.
Zero disables LED-induced tunneling while preserving LED charge resets and
the quiescent contribution. Use `result.led` for pulse and charge-reset truth,
and `result.flip_times` for all recorded parity flips. Each simulation call
starts a fresh chunk with the same LED phase; change the seed for a new
realization. See [`notebooks/led_pulsing.ipynb`](notebooks/led_pulsing.ipynb)
for timing controls, waveform plots, and an optional count sweep.

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

## LED-pulsed records: network reconstruction

For LED-pulsed records, `reconstruct_parity_flips_led` runs a trained network. It has 9848 parameters, ships with
the package and runs on a CPU. It sees the whole record together with the LED schedule,
per-sample summaries of the LED-window HMMs (`qpd.mlebench.led_hmm`) and a record-level burst-intensity estimate,
and returns every flip, both in the quiet periods and inside the bursts. It needs the optional extra:

```bash
pip install -e ".[led]"      # torch, pandas, numba
```

```python
from qpd.reconstruction import reconstruct_parity_flips_led, plot_trace_with_flips, score_flips

rec = reconstruct_parity_flips_led(iq, 1e5, first_pulse_time_s=0.0731, fallback_hmmsf=True)
rec.flip_times    # [s], sorted, after a fixed 300 us veto at each LED onset (the charge reset steps I/Q)
rec.confidence    # per-flip presence probability in [0.5, 1]
rec.degenerate    # True if a hard health flag fired (see rec.flags); fallback_hmmsf then returned HMM-SF instead
rec.mu_hat, rec.rate_hz, rec.hmmsf_flip_times, rec.branch, rec.diagnostics

plot_trace_with_flips(iq, 1e5, rec.flip_times, confidence=rec.confidence, window=(0.072, 0.086))
score_flips(truth, rec.flip_times)

# no trigger time logged: the LED phase is estimated from the waveform (qpd.mlebench.led_phase)
rec = reconstruct_parity_flips_led(iq, 1e5)
rec.first_pulse_time_s, rec.phase_estimated   # the estimate used for every input and for the veto
```

**Without the trigger time.** The phase estimate is within tens of µs at SNR ≥ 2, and F1 is unchanged. At SNR 1,
F1 drops by about 0.01. At SNR 0.5 the estimate is often off by about 1 ms, and F1 drops by 0.02–0.05. Even with the
estimated phase, the network still beats HMM-SF given the true phase in most settings. Pass the trigger time when it
is logged.

- **Schedule and sampling:** the network was trained on 1–5 s records at 100 kHz with the `LEDPulseModel` schedule (50 µs on, 300 ms off), per-sample SNR 0.5–5 and 0–10 tunnelling events per pulse.
- **Where it beats HMM-SF:** inside that range it matches or beats HMM-SF, `predict_led_hmm(..., pin_quiet_rate=True, smooth_fit_s=2e-3, min_smoothed_contrast=2.5)`. The gain is largest when bursts are crowded.
- **Above the trained SNR:** HMM-SF is better and the network degrades. That raises the hard flag `snr_above_trained_range`.
- **Health checks:** a set of flags reports when a record is unlike the training data. Hard flags fire on coloured noise left unwhitened, an incoherent parity axis, no usable LED window, or disagreement with HMM-SF at high SNR.
- **Coloured noise:** it is whitened automatically (`whiten='auto'`).
- **Full example:** [`examples/example_led_network.py`](examples/example_led_network.py); `checks/check_led_network.py` is the gate.

## Verification

```bash
python checks/check_led_pulsing.py              # LED pulse / count / charge checks
python checks/check_led_network.py              # LED flip network (needs .[led])
python checks/check_parity_reconstruction.py     # reconstruction regression gate
python checks/check_reconstruction_benchmark.py  # surrogate-replay benchmark gate
python checks/check_readout_window.py            # dressed-S21 + lock-in window
python checks/study_parity_reconstruction.py     # rate / noise / burst / device sweeps
python checks/study_reconstruction_diagnostics.py  # the three figures above
```

## Notebooks

`notebooks/simulation.ipynb` runs the simulator end to end and, in §7,
demonstrates the reconstruction against the simulator's own ground truth.
