"""Reconstruct parity flips in an LED-pulsed record with the trained network.

Simulates one 2 s LED-pulsed record, runs ``reconstruct_parity_flips_led`` next to the HMM-SF baseline, scores both
against the simulator's truth and plots the first LED pulse.

Needs the optional extra: ``pip install -e ".[led]"``.
Run: ``python examples/example_led_network.py``
"""
import numpy as np
import matplotlib.pyplot as plt

from qpd import QPD, LEDPulseModel, ResonatorConfig, VNASimulator, WhiteGaussianNoise
from qpd.reconstruction import plot_trace_with_flips, score_flips
from qpd.reconstruction.led_network import reconstruct_parity_flips_led

device = QPD(e_j_hz=6.95e9, e_c_hz=0.695e9)
device.coupling_g_hz = 40e6
first_pulse = 0.0731                                  # LED trigger: time of the first pulse onset [s]
sim = VNASimulator(
    qpd=device, resonator=ResonatorConfig(f_r=6.936e9, q_i=1e5, q_c_abs=1.2e5),
    f_drive=6.936437e9, sample_rate=1e5, noise=WhiteGaussianNoise(sigma=0.15),   # per-sample SNR ~2 (trained: 0.5-5)
    led_pulsing=LEDPulseModel(first_pulse_time=first_pulse, expected_n_qp=10),
)
result = sim.simulate(duration=2.0, seed=3)

# fallback_hmmsf: records the network should not decode (rec.degenerate, e.g. SNR above the trained range) get HMM-SF
rec = reconstruct_parity_flips_led(result.iq, 1e5, first_pulse_time_s=first_pulse, fallback_hmmsf=True)
print(f"{rec.flip_times.size} flips  decoder={rec.decoder}  degenerate={rec.degenerate}  flags={rec.flags}")
print(f"mu_hat={rec.mu_hat:.2f} flips per pulse  quiet rate ~ {rec.rate_hz:.1f} Hz")

net = score_flips(result.flip_times, rec.flip_times, tol=200e-6)
sf = score_flips(result.flip_times, rec.hmmsf_flip_times, tol=200e-6)
print(f"network: efficiency {net.efficiency:.3f}  purity {net.purity:.3f}  timing rms {1e6*net.rms_s:.0f} us")
print(f"HMM-SF:  efficiency {sf.efficiency:.3f}  purity {sf.purity:.3f}  timing rms {1e6*sf.rms_s:.0f} us")

fig, _ = plot_trace_with_flips(result.iq, 1e5, rec.flip_times, truth_times=result.flip_times,
                               confidence=rec.confidence, window=(first_pulse - 1e-3, first_pulse + 12e-3),
                               smooth_hz=5e3, title="LED network: first pulse")
plt.show()
