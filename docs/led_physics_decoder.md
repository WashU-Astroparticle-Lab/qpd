# Calibrated LED parity-flip reconstruction

Use this decoder when a readout calibration provides paired even and odd complex I/Q means on an equally weighted offset-charge grid. Both rows must refer to the same charges in the same order. Supply the single-sample noise standard deviation per quadrature in the same I/Q units.

```python
import numpy as np
from qpd.reconstruction import reconstruct_parity_flips_led_physics

result = reconstruct_parity_flips_led_physics(
    I + 1j*Q, sample_rate=10_000, first_pulse_time_s=0.010,
    on_time_s=50e-6, off_time_s=30e-3,
    branch_iq=np.stack([even_iq, odd_iq]), noise_sigma=sigma_iq,
    tolerance_s=200e-6, call_budget=len(native_flip_times),
    return_charge=True,
)
flip_times = result.flip_times         # seconds from sample zero
charge_mass = result.charge_posterior  # [sample, charge bin]
```

The example matches the native call volume for a research comparison. For independent operation supply a positive `min_expected_count` instead of `call_budget`. Calibrate that threshold on independent active and zero-LED records. It is an expected number of flips in a timing window and can exceed one. No validated deployment threshold is provided.

The hidden state comprises parity and offset charge. Each LED rise redraws charge. Additional redraws follow a Poisson process. The flip rate combines quiescent flips with a Poisson yield at each pulse, distributed in time by uniform excitation during ON and a nonnegative-conditioned exponentially modified Gaussian delay. Defaults correspond to one mean LED flip per pulse, a quiet rate of 20 Hz, and a charge redraw rate of 100 per minute. These are model assumptions. The estimator uses the supplied schedule and does not know the realized charge or event history.

Forward–backward inference averages over both charge and parity. `expected_counts[k]` contains eight subinterval expectations between samples k and k+1, including even-multiplicity flips invisible in the sampled parity. Timestamp selection greedily covers disjoint windows of width approximately twice the requested tolerance. Close flips can compete for a window. The optional charge posterior can have multiple modes. With `align_iq=True` the decoder estimates one rigid rotation and translation per record. Input already in the calibration frame needs no alignment.

The likelihood assumes instantaneous parity sampling and independent Gaussian I/Q noise with equal variance in each quadrature. Sample averaging, coloured noise, drift during a record or calibration errors can change the operating point on an instrument.

Working memory grows with samples times charge bins. Run the focused check with `python checks/check_led_physics.py`.
