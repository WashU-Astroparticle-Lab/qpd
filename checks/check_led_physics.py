"""Regression checks for calibrated LED smoothing. Run from the repository root."""
from itertools import product
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from qpd.reconstruction import reconstruct_parity_flips_led_physics as decode


def check_posterior():
    # Enumerate joint states on a short trace with an off-grid scheduled redraw
    # and unscheduled redraws. This checks both backward smoothing and gap counts.
    branches = np.array([[0+0j, 1j], [1+1j, 1+0j]])
    iq = np.array([.1+.1j, .6+.9j, .95+.8j, .2+.8j])
    sigma, fs, rate = .4, 1000, 70
    result = decode(iq, fs, .0013, on_time_s=.00005, off_time_s=.03,
                    branch_iq=branches, noise_sigma=sigma, mean_led_flips=0,
                    quiescent_rate_hz=rate, charge_jump_rate_hz=40,
                    call_budget=0, return_charge=True)
    toggle = -.5*np.expm1(-2*rate/fs)
    charge_sum = np.zeros((len(iq), 2))
    odd_sum = np.zeros(len(iq)-1)
    total = 0.
    for trajectory in product(range(4), repeat=len(iq)):
        parity, charge = np.divmod(trajectory, 2)
        residual = iq - branches[parity, charge]
        weight = np.exp(-.5*(abs(residual)/sigma)**2).prod()/4
        changes = np.diff(parity) != 0
        weight *= np.where(changes, toggle, 1-toggle).prod()
        reset = np.array([-np.expm1(-40/fs), 1., -np.expm1(-40/fs)])
        weight *= (reset/2 + (1-reset)*(np.diff(charge)==0)).prod()
        odd_sum += weight*changes
        charge_sum[np.arange(len(iq)), charge] += weight
        total += weight
    odds = odd_sum/total
    lam = rate/fs
    np.testing.assert_allclose(result.expected_counts.sum(axis=1),
                               lam*((1-odds)*np.tanh(lam)+odds/np.tanh(lam)))
    np.testing.assert_allclose(result.charge_posterior, charge_sum/total)
    assert result.flip_times.size == 0


def check_pulse_and_calls():
    # Flat, indistinguishable branches leave the event intensity equal to its prior.
    kwargs = dict(sample_rate=10000, first_pulse_time_s=.002031,
                  on_time_s=.000047, off_time_s=.2, noise_sigma=.1,
                  branch_iq=np.zeros((2, 3)), quiescent_rate_hz=0,
                  mean_led_flips=2., return_charge=True)
    trace = np.zeros((400, 2))
    result = decode(trace, call_budget=3, tolerance_s=.0002, **kwargs)
    np.testing.assert_allclose(result.expected_counts.sum(), 2., atol=1e-11)
    assert result.expected_counts.ravel()[:162].sum() == 0
    assert len(result.flip_times) == 3 and np.all(np.diff(result.flip_times) >= .0004-1e-14)
    np.testing.assert_allclose(result.charge_posterior, 1/3)
    kwargs.update(first_pulse_time_s=1.)
    quiet = decode(trace, min_expected_count=.01, align_iq=True, **kwargs)
    assert quiet.flip_times.size == 0 and quiet.expected_counts.sum() == 0


if __name__ == "__main__":
    check_posterior()
    check_pulse_and_calls()
    print("LED physics decoder checks passed.")
