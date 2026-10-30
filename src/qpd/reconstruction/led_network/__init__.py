"""Neural-network flip reconstruction for LED-pulsed records (optional extra: ``pip install "qpd[led]"``).

    from qpd.reconstruction.led_network import reconstruct_parity_flips_led
    rec = reconstruct_parity_flips_led(iq, 1e5, first_pulse_time_s=0.1234)   # or omit it: phase estimated
    rec.flip_times, rec.confidence, rec.degenerate

A trained network (9848 values, bundled as ``led_flip_network.pt``) reads the whole record together with per-sample
summaries of the LED-window HMM decoders (:mod:`qpd.mlebench.led_hmm`), the LED schedule and a record-level burst
intensity, and reports every parity flip. It needs torch and pandas; numba makes the burst-intensity estimate fast.
Same call pattern and result fields as the other ``qpd.reconstruction`` entry points.
"""
try:
    import torch  # noqa: F401
except ImportError as error:  # pragma: no cover - depends on the environment
    raise ImportError('qpd.reconstruction.led_network needs torch: pip install "qpd[led]"') from error

from .reconstruct import LEDReconstructionResult, reconstruct_parity_flips_led  # noqa: F401,E402

__all__ = ['reconstruct_parity_flips_led', 'LEDReconstructionResult']
