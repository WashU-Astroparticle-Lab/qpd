#!/usr/bin/env python
"""Gate for offset-charge event detection in the static reconstruction path.

Run: ``python checks/check_charge_event_static.py``; exits non-zero on any
failure. Each numbered section gates one claim from issue #44:

1. **Null calibration** -- on jump-free physics traces the detector accepts
   nothing, and the largest verification gain it ever considers sits well
   below the shipped threshold.
2. **No-jump invariance** -- with no verified event, ``segment_charge_jumps``
   on/off produce identical flip lists: opting in costs nothing on the
   traces where nothing is found.
3. **Mid-trace jump** -- a detectable jump is found, localized, and the flip
   F1 stays at the jump-free level (the "before" decode is printed for
   contrast).
4. **Parity-blind landing** -- the worst measured failure of the stale-model
   decode (F1 0.30, rate inflated 16 -> 40 Hz at delta = 0.07): the post-jump
   stretch must be declared dead, the pre-jump flips kept, the rate sane,
   and no burst fabricated out of the dead window.
5. **Impact event** -- a jump coincident with an n_qp ~ 20 burst (the
   physically correlated case): the burst must survive segmentation, and the
   detected burst is flagged charge-coincident.
6. **Mirror irreducibility** -- delta = 0.5 maps the blob pair onto itself
   (chi_odd(n_g) = chi_even(n_g + 1/2)); the detector must stay silent and
   the decode must be unharmed, because the trace is genuinely unchanged.
7. **API and plumbing** -- diagnostics keys, the result field, benchmark
   passthrough, the ``live=`` mask regression, and the coincidence flag.
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from qpd.reconstruction import (characterize_trace, detect_bursts,  # noqa: E402
                                detect_charge_events,
                                flag_charge_coincidences,
                                reconstruct_parity_flips_static, score_flips)
from qpd.reconstruction.benchmark import _RECON_KEYS  # noqa: E402
from qpd.reconstruction.burst_hmm import decode_burst_aware  # noqa: E402
from qpd.reconstruction.static_bias import fit_two_blobs  # noqa: E402
from reconstruction_scenarios import build_static_scenario  # noqa: E402

NG0 = 0.18
FS = 1e4
BG = 17.0
DUR = 5.0
TOL = 0.5e-3
_failures = []


def check(name, ok, detail=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f"  --  {detail}"
                                                      if detail else ""))
    if not ok:
        _failures.append(name)


def scenario(**kw):
    kw.setdefault("n_g", NG0)
    kw.setdefault("duration", DUR)
    kw.setdefault("sample_rate", FS)
    kw.setdefault("tunnel_rate_hz", BG)
    return build_static_scenario(**kw)


def jump_scenario(delta, t_jump=2.5, **kw):
    return scenario(charge_jumps=(np.array([t_jump]), np.array([delta])), **kw)


# --- 1. null calibration --------------------------------------------------
print("1. Detector null on jump-free traces")
scn = scenario()
n_bound, gains = 0, []
for seed in range(12):
    tr = scn.simulate(seed=7000 + seed)
    d = detect_charge_events(tr.iq, fit_two_blobs(tr.iq), 1.0 / FS)
    n_bound += d.boundaries.size
    gains.append(d.scan["max_gain_nats"])
check("no boundaries accepted on 12 jump-free traces", n_bound == 0,
      f"{n_bound} accepted")
check("largest null verification gain clears the threshold with margin",
      max(gains) < 0.67 * 15.0,
      f"max {max(gains):.1f} vs threshold 15")

# --- 2. no-jump invariance ------------------------------------------------
print("\n2. No-jump invariance of the default")
same = True
for seed in (3, 11):
    tr = scn.simulate(seed=seed)
    r_on = reconstruct_parity_flips_static(tr.iq, FS, segment_charge_jumps=True)
    r_off = reconstruct_parity_flips_static(tr.iq, FS,
                                            segment_charge_jumps=False)
    same &= (r_on.flip_times.size == r_off.flip_times.size
             and np.allclose(r_on.flip_times, r_off.flip_times)
             and r_on.charge_jump_times.size == 0)
check("segment_charge_jumps=True is a no-op when nothing is detected", same)

# --- 3. mid-trace jump ----------------------------------------------------
print("\n3. Mid-trace jump (delta = 0.20)")
scn_j = jump_scenario(0.20)
f1_on, f1_off, errs, found = [], [], [], 0
for seed in range(4):
    tr = scn_j.simulate(seed=seed)
    r_on = reconstruct_parity_flips_static(tr.iq, FS, segment_charge_jumps=True)
    r_off = reconstruct_parity_flips_static(tr.iq, FS,
                                            segment_charge_jumps=False)
    f1_on.append(score_flips(tr.flip_times, r_on.flip_times, tol=TOL).hard_f1)
    f1_off.append(score_flips(tr.flip_times, r_off.flip_times,
                              tol=TOL).hard_f1)
    if r_on.charge_jump_times.size:
        found += 1
        errs.append(float(np.min(np.abs(r_on.charge_jump_times - 2.5))))
check("jump found in every trace", found == 4, f"{found}/4")
check("boundary within 20 ms of truth", max(errs) < 0.02 if errs else False,
      f"worst {1e3 * max(errs):.1f} ms" if errs else "not found")
# Jump-free baseline at this operating point: F1 ~0.95.
check("flip F1 at the jump-free level",
      np.mean(f1_on) > 0.9,
      f"after {np.mean(f1_on):.3f} (before-decode {np.mean(f1_off):.3f})")

# --- 4. parity-blind landing ----------------------------------------------
print("\n4. Blind landing (delta = 0.07): dead time, not silent garbage")
scn_b = jump_scenario(0.07)
ok_dead, ok_pre, rates, fake, rates_off = True, True, [], 0, []
for seed in range(4):
    tr = scn_b.simulate(seed=seed)
    r = reconstruct_parity_flips_static(tr.iq, FS, segment_charge_jumps=True)
    r_off = reconstruct_parity_flips_static(tr.iq, FS,
                                            segment_charge_jumps=False)
    d = r.diagnostics
    dead = d.get("dead_windows", [])
    # the post-jump stretch (2.5 s .. 5 s) must be mostly dead
    dead_span = sum(hi - lo for lo, hi in dead if lo > 2.0)
    ok_dead &= dead_span > 1.5
    pre_truth = tr.flip_times[tr.flip_times < 2.4]
    s_pre = score_flips(pre_truth, r.flip_times[r.flip_times < 2.4], tol=TOL)
    ok_pre &= s_pre.hard_f1 > 0.85
    rates.append(d["rate_hz"])
    rates_off.append(r_off.diagnostics["rate_hz"])
    for b in detect_bursts(r.flip_times, d["rate_hz"], duration=DUR):
        if b.t_start > 2.5:
            fake += 1
check("post-jump stretch declared dead", ok_dead)
check("pre-jump flips still recovered (F1 > 0.85)", ok_pre)
check("rate not inflated by the blind half",
      np.mean(rates) < 25.0,
      f"after {np.mean(rates):.1f} Hz vs before-decode "
      f"{np.mean(rates_off):.1f} Hz (true 17)")
check("no bursts fabricated in the dead window", fake == 0, f"{fake} found")

# --- 5. impact event: jump + burst together -------------------------------
# Two regimes, gated separately. Landing at a *usable* bias (delta = 0.20,
# post-jump contrast ~2.1): the burst must survive segmentation -- the
# boundary guard and the refit must not eat the correlated signal. Landing
# below the usability threshold (delta = 0.10, post-jump contrast ~1.2 --
# the same bias a whole-trace measurement would have flagged degenerate):
# the multiplicity is genuinely unrecoverable there, and the requirement is
# that the loss is *declared* (dead time + the charge event timestamp still
# marks the impact) instead of the pre-upgrade silent undercount.
print("\n5. Impact event (jump + n_qp ~ 20 burst, coincident)")
scn_nb = scenario(burst_times=np.array([2.5]), burst_n_qp=20.0)
scn_jb = jump_scenario(0.20, burst_times=np.array([2.5]), burst_n_qp=20.0)
n_plain, n_jumped, coincident, trials = [], [], 0, 0
for seed in range(6):
    tr0 = scn_nb.simulate(seed=100 + seed)
    tr1 = scn_jb.simulate(seed=100 + seed)
    if not tr0.bursts or not tr1.bursts:
        continue
    trials += 1
    for tr, out in ((tr0, n_plain), (tr1, n_jumped)):
        b = tr.bursts[0]
        r = reconstruct_parity_flips_static(tr.iq, FS, segment_charge_jumps=True)
        out.append(np.count_nonzero(
            (r.flip_times >= b.t_start - 1e-3)
            & (r.flip_times <= b.t_end + 1e-3)))
        if tr is tr1 and r.charge_jump_times.size:
            bursts = detect_bursts(r.flip_times, r.diagnostics["rate_hz"],
                                   duration=DUR)
            flag_charge_coincidences(bursts, r.charge_jump_times, pad=5e-3)
            coincident += any(b_.charge_coincident for b_ in bursts)
check("burst survives a coincident jump to a usable bias",
      np.mean(n_jumped) > 0.8 * np.mean(n_plain),
      f"in-window {np.mean(n_jumped):.1f} with jump vs "
      f"{np.mean(n_plain):.1f} without")
check("coincident bursts get the charge flag",
      coincident >= max(1, trials // 3), f"{coincident}/{trials} flagged")
scn_low = jump_scenario(0.07, burst_times=np.array([2.5]), burst_n_qp=20.0)
declared = 0
for seed in range(3):
    tr = scn_low.simulate(seed=100 + seed)
    r = reconstruct_parity_flips_static(tr.iq, FS, segment_charge_jumps=True)
    d = r.diagnostics
    declared += (r.charge_jump_times.size > 0
                 and d.get("live_fraction", 1.0) < 0.7)
check("blind landing under a burst is declared dead, event timestamped",
      declared == 3, f"{declared}/3 declared")

# --- 6. mirror irreducibility ---------------------------------------------
print("\n6. Mirror jump (delta = 0.50) is invisible and harmless")
scn_m = jump_scenario(0.50)
n_events, dF1 = 0, []
for seed in range(3):
    tr = scn_m.simulate(seed=seed)
    r = reconstruct_parity_flips_static(tr.iq, FS, segment_charge_jumps=True)
    n_events += r.charge_jump_times.size
    f1 = score_flips(tr.flip_times, r.flip_times, tol=TOL).hard_f1
    trc = scn.simulate(seed=seed)
    rc = reconstruct_parity_flips_static(trc.iq, FS, segment_charge_jumps=True)
    f1c = score_flips(trc.flip_times, rc.flip_times, tol=TOL).hard_f1
    dF1.append(f1c - f1)
check("detector silent (the trace distribution is unchanged)", n_events == 0,
      f"{n_events} events")
check("decode unharmed (F1 within 0.05 of jump-free)",
      max(dF1) < 0.05, f"worst deficit {max(dF1):.3f}")

# --- 7. API and plumbing --------------------------------------------------
print("\n7. API and plumbing")
tr = scn_j.simulate(seed=1)
r = reconstruct_parity_flips_static(tr.iq, FS, segment_charge_jumps=True)
d = r.diagnostics
keys = ("charge_event_times", "charge_event_gain_nats",
        "charge_event_localization_s", "segment_edges", "segment_models",
        "dead_windows", "live_fraction", "n_boundary_suppressed",
        "charge_scan", "global_model")
check("segmented diagnostics carry all keys",
      all(k in d for k in keys),
      ", ".join(k for k in keys if k not in d) or "all present")
check("result.charge_jump_times mirrors the ramped API",
      isinstance(r.charge_jump_times, np.ndarray)
      and r.charge_jump_times.size == d["charge_event_times"].size)
fd = characterize_trace(tr.iq, FS, mode="static",
                        segment_charge_jumps=True)
check("TraceFidelity.charge_jump_times populated in static mode",
      fd.charge_jump_times.size == r.charge_jump_times.size,
      f"{fd.charge_jump_times.size} vs {r.charge_jump_times.size}")
check("new kwargs registered in _RECON_KEYS",
      {"charge_min_gain_nats", "charge_min_segment_samples",
       "boundary_guard_samples"} <= _RECON_KEYS)
check("segment_blocks superseded by detected boundaries",
      reconstruct_parity_flips_static(
          tr.iq, FS, segment_blocks=4, segment_charge_jumps=True
      ).diagnostics.get("segment_blocks_superseded") is True)

# live= mask regression: None must reproduce the unmasked decode exactly,
# and an all-live mask must equal it too.
rng = np.random.default_rng(5)
par = (np.cumsum(rng.random(20000) < 1.7e-3) % 2).astype(float)
x = np.where(par > 0, 1.2, -1.2) + rng.normal(0, 1.0, 20000)
mu_a, mu_b = np.full(x.size, 1.2), np.full(x.size, -1.2)
b0 = decode_burst_aware(x, mu_a, mu_b, 1.0, 1e-4)
b1 = decode_burst_aware(x, mu_a, mu_b, 1.0, 1e-4,
                        live=np.ones(x.size, dtype=bool))
check("decode_burst_aware(live=all) identical to live=None",
      b0.p_quiet == b1.p_quiet
      and np.array_equal(b0.parity_path, b1.parity_path))

print()
if _failures:
    print(f"FAILED: {len(_failures)} check(s): {_failures}")
    sys.exit(1)
print("ALL CHECKS PASSED")
