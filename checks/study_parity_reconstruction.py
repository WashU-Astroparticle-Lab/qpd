#!/usr/bin/env python
"""Performance studies for blind parity-flip timing reconstruction.

Run: ``python checks/study_parity_reconstruction.py [study ...]``
with ``study`` one or more of ``rate``, ``burst``, ``device`` (default: all).

  rate    -- efficiency, purity and timing vs the background tunnelling rate.
  burst   -- reconstruction efficiency and timing bias as a function of how
             many quasiparticle tunnels land in the burst window. This is the
             regime where events crowd: two flips inside a few samples are not
             resolvable at all, which sets a ceiling no algorithm can beat, so
             the intrinsic unresolvable fraction is reported alongside.
  device  -- robustness over E_J/E_C and readout noise. Charge dispersion falls
             exponentially with E_J/E_C, so the parity contrast -- and with it
             everything downstream -- shrinks along that axis.

All runs are blind: the reconstruction is given only the I/Q trace and the
sample rate.
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from reconstruction_scenarios import build_scenario  # noqa: E402

from qpd.reconstruction import (burst_report, reconstruct_parity_flips,  # noqa: E402
                               score_flips)

DT = 1e-5


def _run(scn, seed=0):
    result = scn.simulate(seed=seed)
    rec = reconstruct_parity_flips(result.iq, scn.sim.sample_rate)
    return result, rec


def _unresolvable(times, min_gap):
    """Truth events whose nearest neighbour is closer than ``min_gap``."""
    if times.size < 2:
        return 0
    gaps = np.diff(np.sort(times))
    return int(np.count_nonzero(gaps < min_gap))


def study_rate(duration=20.0, seeds=(0, 7)):
    print("=" * 100)
    print("STUDY 1  --  background tunnelling rate (500 Hz ramp, no bursts, "
          "no charge jumps)")
    print("=" * 100)
    print(f"{'rate/Hz':>8s} {'est/Hz':>8s} {'truth':>6s} {'pred':>6s} "
          f"{'effic':>7s} {'purity':>7s} {'hardF1':>7s} {'softF1':>7s} "
          f"{'bias/us':>8s} {'rms/us':>7s} {'contrast':>9s}")
    for rate in (1.0, 3.0, 10.0, 30.0, 100.0, 300.0):
        rows = []
        for seed in seeds:
            scn = build_scenario(duration=duration, tunnel_rate_hz=rate,
                                 burst_rate_hz=0.0, seed=seed)
            res, rec = _run(scn, seed=seed)
            s = score_flips(res.flip_times, rec.flip_times)
            rows.append((s, rec))
        s = rows[0][0]
        eff = np.mean([r[0].efficiency for r in rows])
        pur = np.mean([r[0].purity for r in rows])
        hf1 = np.mean([r[0].hard_f1 for r in rows])
        sf1 = np.mean([r[0].soft_f1 for r in rows])
        bias = np.mean([r[0].bias_s for r in rows]) * 1e6
        rms = np.mean([r[0].rms_s for r in rows]) * 1e6
        est = np.mean([r[1].rate_hz for r in rows])
        con = np.median([np.median(r[1].emission.contrast) for r in rows])
        n_t = int(np.sum([r[0].n_truth for r in rows]))
        n_p = int(np.sum([r[0].n_pred for r in rows]))
        print(f"{rate:8.0f} {est:8.1f} {n_t:6d} {n_p:6d} {eff:7.3f} {pur:7.3f} "
              f"{hf1:7.3f} {sf1:7.3f} {bias:+8.1f} {rms:7.1f} {con:9.2f}")
    print("\nSummed over seeds " + str(seeds) + f"; {duration:g} s per trace.")


def study_burst(duration=20.0, seed=0):
    print()
    print("=" * 100)
    print("STUDY 2  --  efficiency and timing bias vs tunnels in the burst "
          "window")
    print("=" * 100)
    n_qp_grid = np.array([1, 2, 3, 5, 8, 12, 20, 30, 50], dtype=float)
    scn = build_scenario(duration=duration, tunnel_rate_hz=10.0,
                         burst_rate_hz=3.0, burst_n_qp=n_qp_grid, seed=seed)
    res, rec = _run(scn, seed=seed)
    overall = score_flips(res.flip_times, rec.flip_times)
    print(f"whole trace: {overall.n_truth} true flips, "
          f"{overall.n_pred} predicted, efficiency {overall.efficiency:.3f}, "
          f"purity {overall.purity:.3f}, hard F1 {overall.hard_f1:.3f}")
    print(f"{len(res.bursts)} bursts; background tunnelling 10 Hz\n")

    rows = burst_report(res.flip_times, rec.flip_times, res.bursts)
    # Bin by the number of true flips actually inside the window -- that, not
    # the Poisson mean, is what limits what can be recovered.
    bins = [(1, 2), (2, 4), (4, 7), (7, 11), (11, 17), (17, 26), (26, 41),
            (41, 200)]
    print(f"{'flips in window':>16s} {'bursts':>7s} {'truth':>6s} "
          f"{'matched':>8s} {'effic':>7s} {'purity':>7s} {'bias/us':>8s} "
          f"{'rms/us':>7s} {'unresolv':>9s} {'rate/kHz':>9s}")
    for lo, hi in bins:
        sel = [r for r in rows if lo <= r["score"].n_truth < hi]
        if not sel:
            continue
        n_t = sum(r["score"].n_truth for r in sel)
        n_m = sum(r["score"].n_matched for r in sel)
        n_p = sum(r["score"].n_pred for r in sel)
        dts = np.concatenate([r["score"].dt for r in sel
                              if r["score"].dt.size]) if n_m else np.zeros(1)
        unres = sum(_unresolvable(
            res.flip_times[(res.flip_times >= r["t_start"])
                           & (res.flip_times < r["t_end"])], 2 * DT)
            for r in sel)
        rate_khz = np.mean([r["flip_rate_hz"] for r in sel]) / 1e3
        print(f"{f'{lo}-{hi - 1}':>16s} {len(sel):7d} {n_t:6d} {n_m:8d} "
              f"{n_m / max(n_t, 1):7.3f} {n_m / max(n_p, 1):7.3f} "
              f"{np.mean(dts) * 1e6:+8.1f} "
              f"{np.sqrt(np.mean(dts ** 2)) * 1e6:7.1f} "
              f"{unres / max(n_t, 1):9.3f} {rate_khz:9.2f}")
    print("\n'unresolv' = fraction of true flips whose nearest neighbour is "
          "within 2 samples (20 us):\nan intrinsic ceiling on efficiency, not "
          "an algorithm failure. 'rate' is the mean in-window flip rate.")


def study_device(duration=5.0, seed=0):
    print()
    print("=" * 100)
    print("STUDY 3  --  device ratio and readout noise (with bursts and "
          "charge jumps)")
    print("=" * 100)
    print(f"{'E_J/E_C':>8s} {'sigma':>8s} {'contrast':>9s} {'truth':>6s} "
          f"{'pred':>6s} {'effic':>7s} {'purity':>7s} {'hardF1':>7s} "
          f"{'rms/us':>7s} {'comb':>6s}")
    e_c = 0.695e9
    for ratio in (10.0, 12.0, 15.0, 20.0):
        for sigma in (1e-4, 2e-4, 5e-4):
            scn = build_scenario(duration=duration, tunnel_rate_hz=10.0,
                                 e_j_hz=ratio * e_c, e_c_hz=e_c, sigma=sigma,
                                 burst_rate_hz=1.0, charge_jump_rate_hz=0.1,
                                 seed=seed)
            res, rec = _run(scn, seed=seed)
            s = score_flips(res.flip_times, rec.flip_times)
            con = float(np.median(rec.emission.contrast))
            comb = "yes" if rec.reset_comb is not None else "no"
            print(f"{ratio:8.0f} {sigma:8.0e} {con:9.2f} {s.n_truth:6d} "
                  f"{s.n_pred:6d} {s.efficiency:7.3f} {s.purity:7.3f} "
                  f"{s.hard_f1:7.3f} {s.rms_s * 1e6:7.1f} {comb:>6s}")
    print("\n'contrast' = median branch separation in units of the noise sigma.")


def main(argv):
    wanted = [a for a in argv[1:] if not a.startswith("-")] or \
        ["rate", "burst", "device"]
    if "rate" in wanted:
        study_rate()
    if "burst" in wanted:
        study_burst()
    if "device" in wanted:
        study_device()
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
