"""Adapter example. Native HMM-SF here gives a reproducible equality check."""
from qpd.mlebench.led_generalization_hmmsf import predict_hmmsf


def reconstruct(iq, *, sample_rate_hz, first_pulse_time_s, on_time_s, off_time_s):
    result = predict_hmmsf(
        iq, sample_rate_hz, first_pulse_time_s=first_pulse_time_s,
        on_time_s=on_time_s, off_time_s=off_time_s,
        nominal_emg={"mu_s": .0002, "sigma_s": .0004, "tau_s": .0008},
    )
    return result.flip_times
