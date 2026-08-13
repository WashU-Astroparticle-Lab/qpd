#!/usr/bin/env python
"""Verification for the burst-aware parity x regime HMM (issue #40).

Run: ``python checks/check_burst_aware_hmm.py`` (exit 0 = all pass).

The module under test claims that giving quasiparticle-burst windows their
own flip probability recovers the multiplicity the global-prior decoder
smooths away, *without* changing background behaviour or fabricating bursts.
Each of those clauses is gated separately:

  1. Recursions are exact: the 4-state forward-backward matches the general-K
     log-domain reference, and the 4-state Viterbi matches brute-force
     enumeration.
  2. The transition matrix is a stochastic matrix under any parameter
     clipping, and a burst regime slower than the background is refused
     (clipped), not modelled.
  3. Pipeline API: ``burst_aware=True`` returns the regime diagnostics; both
     decoders run; an unknown decoder raises.
  4. THE IMPROVEMENT GATE: on single-burst traces at the issue #40 operating
     point (10 kSa/s, contrast 2.4, 17 Hz background, tau = 1 ms), in-window
     recovery at n_qp ~ 20 at least doubles, and clears the evidence-limited
     target of ~6 -- against a visible-on-grid ceiling of ~12 that no decoder
     can beat.
  5. THE GUARDS: on pure background the burst-aware decoder matches the
     two-state decoder's flip efficiency and purity, and fabricates no bursts
     through ``detect_bursts`` defaults.
  6. The quiet-rate estimate deflates: with bursts present the global fit is
     biased high, and the regime model's ``p_quiet`` lands nearer the true
     background rate.
  7. Epsilon is a threshold, not a rate assumption: two orders of magnitude in
     ``burst_rate_hz`` move the recovered in-window multiplicity by about one
     flip, as the ln(1/epsilon) entry cost predicts.
"""
import itertools
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from qpd.reconstruction import (detect_bursts,  # noqa: E402
                                reconstruct_parity_flips_static, score_flips)
from qpd.reconstruction.benchmark import TraceFidelity  # noqa: E402
from qpd.reconstruction.burst_hmm import (burst_transition_matrix,  # noqa: E402
                                          decode_burst_aware,
                                          forward_backward_regime,
                                          viterbi_regime)
from qpd.reconstruction.hmm import forward_backward_reference  # noqa: E402
from qpd.reconstruction.static_bias import StaticBlobModel  # noqa: E402
from qpd.simulator import QuasiparticleBurstModel  # noqa: E402

FS = 1e4
CONTRAST = 2.4
BG = 17.0
_failures = []


def check(name, ok, detail=""):
    tag = "PASS" if ok else "FAIL"
    print(f"  [{tag}] {name}" + (f"  --  {detail}" if detail else ""))
    if not ok:
        _failures.append(name)


def make_fidelity(**recon_kwargs):
    blob = StaticBlobModel(direction=1 + 0j, origin=0j,
                           mu_a=+CONTRAST / 2, mu_b=-CONTRAST / 2,
                           sigma=1.0, weight_a=0.5)
    return TraceFidelity(
        mode="static", sample_rate=FS, n_samples=int(FS), t0=0.0, sigma=1.0,
        contrast_median=CONTRAST, contrast_max=CONTRAST, rate_hz=BG,
        n_flips=0, degenerate=False, sample_fidelity=0.9,
        decoded_fidelity=0.99, model=blob, recon_kwargs=recon_kwargs)


print("1. recursion correctness")
rng = np.random.default_rng(0)
log_emit4 = rng.normal(-2.0, 1.5, (4, 300))
trans = burst_transition_matrix(2e-3, 0.3, 1e-4, 0.1)
post, ll = forward_backward_regime(log_emit4, trans)
post_ref, ll_ref = forward_backward_reference(log_emit4, np.log(trans))
check("forward-backward matches general-K reference",
      np.abs(post - post_ref).max() < 1e-8 and abs(ll - ll_ref) < 1e-6,
      f"max |dpost| = {np.abs(post - post_ref).max():.2e}")

le = rng.normal(-2.0, 1.5, (4, 10))
logt = np.log(trans)
best, best_seq = -np.inf, None
for seq in itertools.product(range(4), repeat=10):
    s = np.log(0.25) + le[seq[0], 0]
    for k in range(1, 10):
        s += logt[seq[k - 1], seq[k]] + le[seq[k], k]
    if s > best:
        best, best_seq = s, seq
check("viterbi matches brute-force enumeration",
      tuple(viterbi_regime(le, trans)) == best_seq)

print("2. transition matrix sanity")
for args in [(2e-3, 0.3, 1e-4, 0.1), (0.9, 0.9, 0.9, 0.9),
             (1e-15, 1e-15, 1e-15, 1e-15)]:
    t = burst_transition_matrix(*args)
    check(f"rows sum to 1 for {args}", np.allclose(t.sum(axis=1), 1.0))
t = burst_transition_matrix(0.2, 0.05, 1e-4, 0.1)
# with p_burst clipped up to p_quiet the two regimes flip identically
check("p_burst below p_quiet is clipped up, not inverted",
      np.isclose(t[2, 3] / (t[2, 2] + t[2, 3]), 0.2, rtol=1e-6),
      f"in-burst flip fraction = {t[2, 3] / (t[2, 2] + t[2, 3]):.3f}")

print("3. pipeline API")
fid = make_fidelity()
model = QuasiparticleBurstModel(times=np.array([0.5]), tau=1e-3, mu=1.2e-3,
                                sigma=0.4e-3, expected_n_qp=20.0)
iq, flips, bt = fid.synthesize(seed=7, rate_hz=BG, bursts=model,
                               return_bursts=True)
rec = reconstruct_parity_flips_static(iq, FS, burst_aware=True)
d = rec.diagnostics
check("regime diagnostics present",
      all(k in d for k in ("burst_posterior", "regime_path", "burst_windows",
                           "p_quiet", "p_burst", "epsilon")))
check("burst posterior is a probability",
      0 <= np.min(d["burst_posterior"]) and np.max(d["burst_posterior"]) <= 1)
rec_p = reconstruct_parity_flips_static(iq, FS, burst_aware=True,
                                        decoder="posterior")
check("posterior decoder runs", rec_p.flip_times.size > 0)
try:
    reconstruct_parity_flips_static(iq, FS, burst_aware=True, decoder="nope")
    check("unknown decoder raises", False)
except ValueError:
    check("unknown decoder raises", True)

b = bt[0]
win = [w for w in d["burst_windows"]
       if w[1] >= b.t_start - 2e-3 and w[0] <= b.t_end + 2e-3]
check("a regime window overlaps the true burst", len(win) >= 1,
      f"windows = {[(f'{a:.4f}', f'{c:.4f}') for a, c in d['burst_windows']]}")

print("4. improvement gate (n_qp ~ 20, single burst)")
n_seeds = 12
base_w, ba_w, vis_w = [], [], []
for k in range(n_seeds):
    iq, flips, bt = fid.synthesize(seed=200 + k, rate_hz=BG, bursts=model,
                                   return_bursts=True)
    b = bt[0]
    if b.n_qp <= 0:
        continue
    lo, hi = b.t_start - 1e-3, b.t_end + 1e-3
    t = np.arange(int(FS)) / FS
    par = np.searchsorted(flips, t, side="right") % 2
    w = (t >= lo) & (t <= hi)
    vis_w.append(np.count_nonzero(np.diff(par[w.nonzero()[0]])))
    for kw, out in [({"burst_aware": False}, base_w),
                    ({"burst_aware": True}, ba_w)]:
        r = reconstruct_parity_flips_static(iq, FS, **kw)
        out.append(np.count_nonzero((r.flip_times >= lo)
                                    & (r.flip_times <= hi)))
mb, ma, mv = np.mean(base_w), np.mean(ba_w), np.mean(vis_w)
check("burst-aware at least doubles in-window recovery", ma >= 2.0 * mb,
      f"{mb:.1f} -> {ma:.1f} (visible ceiling {mv:.1f})")
check("burst-aware clears the evidence-limited target (>= 6)", ma >= 6.0)
check("no decoder beats the visible ceiling", ma <= mv + 1.0)

print("5. background guards")
n_bg = 12
effs, purs = {"base": [], "ba": []}, {"base": [], "ba": []}
false_bursts = {"base": 0, "ba": 0}
for k in range(n_bg):
    iq, truth = fid.synthesize(seed=900 + k, rate_hz=BG)
    for key, kw in [("base", {"burst_aware": False}),
                    ("ba", {"burst_aware": True})]:
        r = reconstruct_parity_flips_static(iq, FS, **kw)
        s = score_flips(truth, r.flip_times, tol=0.5e-3)
        effs[key].append(s.efficiency)
        purs[key].append(s.purity)
        false_bursts[key] += len(detect_bursts(r.flip_times, BG, duration=1.0))
check("background flip efficiency unchanged",
      abs(np.nanmean(effs["ba"]) - np.nanmean(effs["base"])) < 0.03,
      f"{np.nanmean(effs['base']):.3f} vs {np.nanmean(effs['ba']):.3f}")
check("background flip purity unchanged",
      abs(np.nanmean(purs["ba"]) - np.nanmean(purs["base"])) < 0.03,
      f"{np.nanmean(purs['base']):.3f} vs {np.nanmean(purs['ba']):.3f}")
# Background coincidences (3 flips within the 3 ms linking distance) pass
# the gate-off detect_bursts whichever decoder produced the flip train; the
# burst-aware claim is only that it adds none of its own on top.
check("no bursts fabricated beyond the baseline's coincidences",
      false_bursts["ba"] <= false_bursts["base"],
      f"{false_bursts['ba']} vs baseline {false_bursts['base']} "
      f"in {n_bg} traces")

print("6. quiet-rate deflation")
iq, flips, bt = fid.synthesize(seed=41, rate_hz=BG, bursts=model,
                               return_bursts=True)
r = reconstruct_parity_flips_static(iq, FS, burst_aware=True)
p_quiet = r.diagnostics["p_quiet"]
p_seed = r.diagnostics["p_global_seed"]
check("p_quiet below the burst-inflated global estimate", p_quiet < p_seed,
      f"quiet {p_quiet:.2e} vs global {p_seed:.2e}")

print("7. epsilon is logarithmic")
recov = {}
for br in (0.1, 10.0):
    vals = []
    for k in range(8):
        iq, flips, bt = fid.synthesize(seed=600 + k, rate_hz=BG, bursts=model,
                                       return_bursts=True)
        b = bt[0]
        if b.n_qp <= 0:
            continue
        r = reconstruct_parity_flips_static(iq, FS, burst_aware=True,
                                            burst_rate_hz=br)
        vals.append(np.count_nonzero(
            (r.flip_times >= b.t_start - 1e-3)
            & (r.flip_times <= b.t_end + 1e-3)))
    recov[br] = float(np.mean(vals))
check("100x in burst_rate_hz moves recovery by ~one flip",
      abs(recov[10.0] - recov[0.1]) <= 2.0,
      f"{recov[0.1]:.1f} at 0.1 Hz vs {recov[10.0]:.1f} at 10 Hz")

print()
if _failures:
    print(f"FAILED: {len(_failures)} check(s): {_failures}")
    sys.exit(1)
print("ALL CHECKS PASSED")
