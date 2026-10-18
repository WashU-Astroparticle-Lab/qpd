#!/usr/bin/env python
"""Verify LED pulse timing, charge resets, count statistics, and VNA integration.

Run: ``python checks/check_led_pulsing.py`` (exit 0 = all checks pass).
Statistical checks use fixed seeds and generous multi-sigma tolerances.
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from qpd import QPD  # noqa: E402
from qpd.simulator import (  # noqa: E402
    ConstantNg,
    LEDPulseModel,
    LEDTruth,
    QuasiparticleBurstModel,
    ResonatorConfig,
    VNASimulator,
    led_count_sweep,
)

_failures = []


def check(name, ok, detail=""):
    tag = "PASS" if ok else "FAIL"
    print(f"  [{tag}] {name}" + (f"  --  {detail}" if detail else ""))
    if not ok:
        _failures.append(name)


def rejects(name, action):
    try:
        action()
    except ValueError:
        check(name, True)
    else:
        check(name, False)


class FastSimulator(VNASimulator):
    """Use inexpensive, charge-dependent shifts for stochastic integration checks."""

    def _chi_grid(self, n_g_lo, n_g_hi):
        grid = np.linspace(n_g_lo - 0.01, n_g_hi + 0.01, 101)
        shift = 2e4 * np.cos(2 * np.pi * grid)
        return grid, shift, -shift


def same_charge(a, b):
    return a.initial_n_g == b.initial_n_g and all(
        np.array_equal(getattr(a, field), getattr(b, field))
        for field in (
            "charge_reset_times", "charge_reset_values", "charge_reset_is_led"
        )
    )


def main():
    print("1. DAQ alignment, zero-count pulses, and absolute charge resets")
    zero_model = LEDPulseModel(expected_n_qp=0)
    zero = zero_model.sample(5.0, np.random.default_rng(121))
    expected_times = 0.05 + np.arange(17) * 0.30005
    check("public truth type", isinstance(zero, LEDTruth))
    check("17 pulses at 50 ms + k * 300.05 ms",
          np.allclose(zero.pulse_times, expected_times, atol=1e-14, rtol=0))
    check("each ON window lasts 50 us",
          np.allclose(zero.pulse_end_times - zero.pulse_times, 50e-6))
    check("zero excitation retains every pulse and reset",
          len(zero.bursts) == 17 and zero.event_times.size == 0
          and all(b.n_qp == 0 for b in zero.bursts)
          and np.array_equal(zero.charge_reset_times[zero.charge_reset_is_led],
                             zero.pulse_times))
    reset_t = zero.charge_reset_times
    before = np.nextafter(reset_t, -np.inf)
    previous_values = np.r_[zero.initial_n_g, zero.charge_reset_values[:-1]]
    check("charge reset is absolute and applies exactly at its timestamp",
          np.array_equal(zero.evaluate_n_g(reset_t), zero.charge_reset_values)
          and np.array_equal(zero.evaluate_n_g(before), previous_values))
    midpoints = np.r_[reset_t[:-1] + np.diff(reset_t) / 2, (reset_t[-1] + 5) / 2]
    check("charge remains fixed between resets",
          np.array_equal(zero.evaluate_n_g(midpoints), zero.charge_reset_values))
    charge = np.r_[zero.initial_n_g, zero.charge_reset_values]
    check("default charge values lie in [-0.5, 0.5)",
          np.all(charge >= -0.5) and np.all(charge < 0.5))
    empty = zero_model.sample(0, np.random.default_rng(121))
    check("zero duration has no pulses or resets",
          empty.pulse_times.size == empty.charge_reset_times.size == 0
          and empty.evaluate_n_g(np.array([0.0]))[0] == empty.initial_n_g)
    at_end = zero_model.sample(0.05, np.random.default_rng(121))
    check("pulse exactly at chunk end is excluded", at_end.pulse_times.size == 0)

    print("2. Independent Poisson charge jumps and uniform redraws")
    jump_counts = []
    redraws = []
    rng = np.random.default_rng(982)
    for _ in range(1800):
        truth = zero_model.sample(5.0, rng)
        spontaneous = ~truth.charge_reset_is_led
        jump_counts.append(np.count_nonzero(spontaneous))
        redraws.extend(truth.charge_reset_values)
    jump_counts = np.asarray(jump_counts)
    total = int(jump_counts.sum())
    expected_total = jump_counts.size * 5 / 60
    check("spontaneous jump rate is one per minute",
          abs(total - expected_total) < 6 * np.sqrt(expected_total),
          f"{total} jumps in {5 * jump_counts.size} s, expected {expected_total:g}")
    check("spontaneous count variance agrees with Poisson mean",
          abs(jump_counts.var() - 5 / 60) < 0.03,
          f"mean={jump_counts.mean():.4f}, variance={jump_counts.var():.4f}")
    redraws = np.asarray(redraws)
    check("charge redraws have uniform mean and variance",
          abs(redraws.mean()) < 0.01 and abs(redraws.var() - 1 / 12) < 0.005,
          f"mean={redraws.mean():.4f}, variance={redraws.var():.4f}")

    print("3. Poisson pulse counts, causal delays, and finite ON windows")
    count_truth = LEDPulseModel(expected_n_qp=15).sample(
        1200, np.random.default_rng(291))
    counts = np.asarray([b.n_qp for b in count_truth.bursts])
    check("per-pulse count mean and variance agree with the configured mean",
          abs(counts.mean() - 15) < 0.4 and abs(counts.var() - 15) < 1.5,
          f"{counts.size} pulses: mean={counts.mean():.3f}, variance={counts.var():.3f}")
    check("all default tunneling events follow their pulse onset",
          all(np.all(b.event_times >= b.t_arrival) for b in count_truth.bursts))
    check("truth preserves each count and exact burst span",
          all(b.n_qp == len(b.event_times) and b.t_start == b.event_times[0]
              and b.t_end == b.event_times[-1]
              for b in count_truth.bursts if b.n_qp))
    uniform_model = LEDPulseModel(
        on_time=0.02, off_time=1, expected_n_qp=10000, tau=0, mu=0, sigma=0)
    uniform = uniform_model.sample(0.1, np.random.default_rng(201)).bursts[0]
    offsets = uniform.event_times - uniform.t_arrival
    check("zero-delay excitation is spread uniformly across the ON window",
          np.all((offsets >= 0) & (offsets < 0.02))
          and abs(offsets.mean() - 0.01) < 0.0004
          and abs(offsets.var() - 0.02**2 / 12) < 0.000003,
          f"delay mean={offsets.mean():.6f} s, variance={offsets.var():.8f} s^2")
    causal_model = LEDPulseModel(
        on_time=1e-12, off_time=1, expected_n_qp=10000,
        tau=0, mu=0, sigma=0.003)
    causal = causal_model.sample(0.1, np.random.default_rng(420)).bursts[0]
    delays = causal.event_times - causal.t_arrival
    # A Gaussian conditioned to be positive has a half-normal mean. Clamping
    # negative delays to zero would instead halve the expected mean.
    check("negative Gaussian delays are conditioned away, not clipped",
          np.all(delays >= 0)
          and abs(delays.mean() - 0.003 * np.sqrt(2 / np.pi)) < 0.0001
          and np.count_nonzero(delays < 1e-8) < 5,
          f"mean delay={delays.mean():.6f} s")

    print("4. Reproducibility and independent charge random streams")
    model = LEDPulseModel(expected_n_qp=31, charge_jump_rate_hz=2)
    a = model.sample(5, np.random.default_rng(71))
    b = model.sample(5, np.random.default_rng(71))
    c = LEDPulseModel(expected_n_qp=0, charge_jump_rate_hz=2).sample(
        5, np.random.default_rng(71))
    d = model.sample(5, np.random.default_rng(72))
    check("same seed reproduces pulse counts, event times, and charge",
          np.array_equal(a.event_times, b.event_times) and same_charge(a, b))
    check("changing LED count does not alter charge resets", same_charge(a, c))
    check("different seeds retain DAQ alignment but redraw chunk state",
          np.array_equal(a.pulse_times, d.pulse_times)
          and not same_charge(a, d) and not np.array_equal(a.event_times, d.event_times))
    sweep = led_count_sweep()
    check("default count sweep is zero plus 31 logarithmic means",
          sweep.size == 32 and sweep[0] == 0 and sweep[1] == 1
          and sweep[-1] == 1000
          and np.allclose(np.diff(np.log10(sweep[1:])), 0.1))

    print("5. VNA integration, quiescent flips, and independent chunks")
    qpd = QPD(e_j_hz=6.95e9, e_c_hz=0.695e9, coupling_g_hz=40e6)
    resonator = ResonatorConfig(f_r=6.936e9, q_i=1e5, q_c_abs=1.2e5)
    base = dict(qpd=qpd, resonator=resonator, f_drive=6.936e9, sample_rate=1000)
    sim = FastSimulator(**base, led_pulsing=zero_model)
    result = sim.simulate(seed=41)
    again = sim.simulate(seed=41)
    check("default duration is five seconds and both quiescent rates are 20 Hz",
          result.t.size == 5000 and sim.gamma_even_to_odd == 20
          and sim.gamma_odd_to_even == 20)
    check("waveform includes LED truth and evaluates its charge path",
          isinstance(result.led, LEDTruth) and result.bursts is result.led.bursts
          and np.array_equal(result.n_g, result.led.evaluate_n_g(result.t)))
    check("repeated simulations are independent chunks with reproducible seeds",
          np.array_equal(result.iq, again.iq)
          and np.array_equal(result.flip_times, again.flip_times)
          and np.array_equal(result.led.pulse_times, expected_times))
    quiet_counts = np.asarray([sim.simulate(seed=k).flip_times.size for k in range(100)])
    check("zero LED count retains the independent 20 Hz parity background",
          abs(quiet_counts.mean() - 20 * result.t[-1]) < 6,
          f"mean={quiet_counts.mean():.2f} flips per 5 s (expected about 100)")

    print("6. Complete final-pulse truth and observable waveform boundaries")
    tail_model = LEDPulseModel(
        first_pulse_time=0.049, on_time=0.01, off_time=1,
        expected_n_qp=1000, tau=0.02, mu=0, sigma=0, charge_jump_rate_hz=0)
    tail_sim = FastSimulator(**base, led_pulsing=tail_model,
                             gamma_even_to_odd=0, gamma_odd_to_even=0)
    tail = tail_sim.simulate(duration=0.055, seed=8)
    truth_events = tail.led.event_times
    visible = truth_events[(truth_events > tail.t[0]) & (truth_events <= tail.t[-1])]
    check("a final partial ON window retains the full pulse count and tail",
          tail.led.pulse_end_times[-1] > 0.055
          and tail.bursts[-1].n_qp == truth_events.size
          and np.any(truth_events > 0.055)
          and np.any(truth_events < tail.t[-1]))
    check("waveform flip truth contains only recorded events",
          np.array_equal(tail.flip_times, visible))
    expected_parity = tail.parity[0] ^ (np.searchsorted(visible, tail.t, side="right") % 2)
    check("waveform parity follows every visible event, including sub-sample flips",
          np.array_equal(tail.parity, expected_parity))

    print("7. Legacy behavior and a real QPD dispersive integration")
    old_default = FastSimulator(**base).simulate(duration=0.1, seed=6)
    old_explicit = FastSimulator(**base, offset_charge=ConstantNg(0.25)).simulate(
        duration=0.1, seed=6)
    check("non-LED mode retains the constant charge and seeded legacy behavior",
          old_default.led is None and old_default.bursts == []
          and np.all(old_default.n_g == 0.25)
          and np.array_equal(old_default.iq, old_explicit.iq)
          and np.array_equal(old_default.flip_times, old_explicit.flip_times))
    real = VNASimulator(**base, led_pulsing=LEDPulseModel(expected_n_qp=4),
                        n_g_grid_points=11, charge_cutoff=8).simulate(
                            duration=0.1, seed=55)
    check("real dispersive calculation produces finite charge-dependent I/Q",
          real.iq.shape == (100,) and np.all(np.isfinite(real.iq))
          and real.led.pulse_times.size == 1 and np.ptp(real.iq.real) > 0)

    print("8. Invalid settings and incompatible simulation modes")
    for kwargs in (
        {"on_time": 0}, {"off_time": -1}, {"first_pulse_time": -1},
        {"expected_n_qp": -1}, {"expected_n_qp": np.nan}, {"tau": -1},
        {"sigma": np.inf}, {"charge_jump_rate_hz": -1},
        {"n_g_min": 1, "n_g_max": 1},
    ):
        rejects(f"rejects {kwargs}", lambda kwargs=kwargs: LEDPulseModel(**kwargs))
    for duration in (-1, np.nan, np.inf):
        rejects(f"rejects duration {duration}",
                lambda duration=duration: zero_model.sample(duration, np.random.default_rng(1)))
    for points in (0, -1, 1.5, True):
        rejects(f"rejects sweep density {points}",
                lambda points=points: led_count_sweep(points))
    rejects("LED mode rejects an explicit offset-charge model",
            lambda: FastSimulator(**base, led_pulsing=zero_model,
                                  offset_charge=ConstantNg(0.25)).simulate(0.1))
    burst = QuasiparticleBurstModel(np.array([0.05]), tau=0.001, mu=0, sigma=0)
    rejects("LED mode rejects a simultaneous standalone burst model",
            lambda: FastSimulator(**base, led_pulsing=zero_model,
                                  quasiparticle_bursts=burst).simulate(0.1))

    print()
    if _failures:
        print(f"RESULT: {len(_failures)} check(s) FAILED: {_failures}")
        return 1
    print("RESULT: all checks PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
