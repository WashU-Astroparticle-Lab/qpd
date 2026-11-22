# Testing a parity-flip reconstruction with QPD

Your algorithm takes one measured complex I/Q record and its acquisition timing and returns flip times in seconds. QPD can generate the records, score every physical flip against your output, and draw the same nominal waveform for every method. Start by running the supplied native HMM-SF adapter as the candidate. The two outputs should agree exactly. Then substitute your algorithm.

| Question | Data to generate | Result to use |
|---|---|---|
| How does my algorithm perform over the campaign’s range of conditions? | Mixed train and validation splits | Six-stratum C1 score on mixed validation |
| Can I train and validate at each exact LED mean and SNR? | Controlled train and validation splits | Independent waveforms at each setting, plus a matched nominal active/zero comparison |
| How do efficiency and purity depend on LED mean or SNR? Where does my algorithm miss flips? | Recurring controlled display split | Two sweep figures and one fixed 0.5-second waveform |

## 1. Run the example

Run these commands from the root of **this PR's checkout**. Use Python 3.11 in the `straxion` conda environment. On a new machine, create it first with `conda create -n straxion python=3.11`. The example runs on a CPU.

```bash
conda activate straxion
python -m pip install -e . 'pandas>=1.5' 'matplotlib>=3.6'

D=led-work
A=examples/led_candidate.py
python scripts/prepare_led_controlled.py --split display --replicates 2 --out "$D/display-smoke"
python scripts/evaluate_led_candidate.py --data "$D/display-smoke" --candidate "$A" \
    --label "Adapter check" --out "$D/smoke-report"
python scripts/plot_led_nominal_diagnostics.py --data "$D/smoke-report" --out "$D/smoke-report/figures"
```

This runs 32 one-second records across 16 conditions. Open `waveform_nominal_0p5s.png` under the figures directory. The first half-second has **21 true flips and 10 native HMM-SF calls**, with 9 strict matches at 200 µs. The adapter and native rows should line up. Two records per point exercise the workflow. Use the 64-record display below for reporting. Give generators and evaluators a fresh output directory each time.

## 2. Generate training and validation data

**For a separate sample at every fixed setting**, choose a replicate count per split. These commands generate 16 settings × 64 training records and 16 × 16 validation records. Set `--cross` on both commands to generate the full 9 × 8 grid (72 settings, including zero LED at each SNR). Without that flag, you get the two one-parameter sweeps with their common nominal point counted once.

```bash
python scripts/prepare_led_controlled.py --split train --replicates 64 --out "$D/fixed-train"
python scripts/prepare_led_controlled.py --split validation --replicates 16 --out "$D/fixed-val"
```

| Variable swept | Exact settings | Other axis in the one-parameter sweep |
|---|---|---|
| Mean LED-induced flips per pulse | 0, 0.1, 0.3, 0.5, 1, 2, 3, 5, 10 | Maximum SNR = 2 |
| Maximum single-sample SNR | 0.5, 0.75, 1, 1.5, 2, 3, 4, 5 | Mean flips per pulse = 1 |

Fixed-condition train and validation use independent seeds per split, setting and replicate. Repeating a command reproduces the same records. Increasing the count extends a split. These panels use a 10 ms first onset and the nominal device response. A zero LED mean disables induced flips. Background parity flips and the charge redraw at each scheduled rise still occur.

**For mixed-condition training and the primary C1 comparison**, use the scientific generator. Omitting `--out` prints the design counts and raw I/Q bytes before you allocate storage. The full design uses 20,000 train and 4,000 validation records. Start with a stratified prefix of 200 records to check your pipeline.

```bash
python scripts/prepare_led_generalization_scientific.py --split train
python scripts/prepare_led_generalization_scientific.py --split train --count 200 --out "$D/mixed-train-200"
python scripts/prepare_led_generalization_scientific.py --split validation --count 200 --out "$D/mixed-val-200"

# Full splits, when ready
python scripts/prepare_led_generalization_scientific.py --split train --out "$D/mixed-train"
python scripts/prepare_led_generalization_scientific.py --split validation --out "$D/mixed-val"
```

Each 200-record block contains 120 broad, 50 nominal and 30 stress records, with exactly 10% zero LED in each component. Broad records sample the ranges below by Latin hypercube. The nominal proposal alternates exact values and clipped ±20% jitter. Stress draws fix combinations of endpoints for low-rate dense pulses, noisy frequent resets, short sparse records, maximum record size, broad overlapping tails, and narrow bursts with high background. Other parameters retain their broad draws. The first onset varies uniformly over one LED period. A random global phase and complex offset are applied to the complete noisy waveform.

| Parameter | Nominal | Mixed design range |
|---|---:|---:|
| Sampling rate | 10 kHz | 5–500 kHz |
| Duration | 1 s | 0.5–5 s |
| LED ON / OFF | 50 µs / 30 ms | 1–500 µs / 5–100 ms |
| Mean LED-induced flips per pulse | 1 | 0–10 |
| Maximum SNR | 2 | 0.5–5 |
| Tunneling delay μ / σ / τ | 0.2 / 0.4 / 0.8 ms | μ fixed, σ 0.2–0.6 ms, τ 0.4–1.2 ms |
| Background parity-flip rate | 20 Hz | 1–100 Hz |
| Additional charge redraws | 100/min | 60–600/min |

Maximum SNR is the largest noiseless parity separation over offset charge divided by the white-noise standard deviation in one quadrature of one native sample. The realized contrast varies with charge. The delay parameters describe the tunneling times after excitation during LED ON. Mixed train and validation both exclude two joint holdouts, (rate ≤10 kHz, SNR ≤1.4, τ ≥0.85 ms) and (rate ≥100 kHz, OFF ≤10 ms, τ ≥1 ms). The protected reserve stays outside this workflow. Exact parameter values and sampling rules are in the [active SNR 2 specification](../configs/siderius/led_standalone_snr2.json).

**Load one labeled training example** with the same reader for either generator. Each directory has a public manifest with opaque IDs and paths relative to `public/`. Per-record files store float32 I and Q. The exact physical event times, including unresolved multiplicities and background events, are in `flip_times`.

```python
import json
from pathlib import Path
import numpy as np

split = Path("led-work/fixed-train")
rows = json.loads((split / "public/manifest.json").read_text())
row = rows[0]
with np.load(split / "public" / row["file"]) as waveform:
    iq = waveform["I"] + 1j * waveform["Q"]
with np.load(split / "private/truth" / (row["record_id"] + ".npz")) as labels:
    target_seconds = labels["flip_times"]
```

To train a separate model at one setting, select training IDs from `private/design.json` by `parameters.mean_tunnels` and `parameters.snr_max`, then join to the public manifest by `record_id`. Use the same filter on the independent validation split. A single blind reconstruction gets only I/Q, sample rate, first onset and LED ON/OFF durations. Keep the condition, truth, seeds, stratum and charge information on the training-label or scoring side. Split complete records before deriving crops. Freeze the checkpoint and thresholds on a fitting/calibration partition before scoring the evaluation records.

## 3. Score your algorithm against native HMM-SF

Copy the [example adapter](../examples/led_candidate.py) into your own Python file and replace the body of `reconstruct(iq, *, sample_rate_hz, first_pulse_time_s, on_time_s, off_time_s)`. Load trained weights once at module import. Return a one-dimensional sequence of finite timestamps in `[0, len(iq)/sample_rate_hz)`, in seconds from the first sample. Return an empty sequence for no detections. The driver sorts timestamps and preserves duplicates. It calls each method on the same complete record.

```bash
A=examples/led_candidate.py   # Point this at your adapter
python scripts/evaluate_led_candidate.py --candidate "$A" --data "$D/mixed-val-200" --out "$D/c1"
python scripts/evaluate_led_candidate.py --candidate "$A" --data "$D/fixed-val" --out "$D/fixed-scores"
```

Read `c1_candidate.json` beside `c1_hmmsf.json`. **Lower C1 error is better.** Whole-record one-to-one assignment minimizes timing cost + 3 per unmatched prediction + 1 per unmatched truth. Timing cost is zero inside ε = max(50 µs, 2/sample rate), then rises quadratically to 4 at the matching cutoff, max(2 ms, 4ε). Pairs at the cutoff remain unmatched. There is no onset veto. Native HMM-SF uses the public pulse schedule at the native sample rate, 2 ms smoothed emission fits, a 2.5 contrast gate and the shared nominal τ.

The primary `error` pools total cost divided by total truth count within each of six strata, then weights broad/nominal/stress by 0.60/0.25/0.15 and zero/active by 0.10/0.90. Check `primary_valid` first. A missing or wholly zero-truth stratum gives `error: null`, as expected on fixed-condition panels. Use `by_record`, `by_stratum`, `pooled` and `empty_truth` for precise and late matches, fixed-50-µs matches, FP, FN, cost/s and FP/s. Empty records still incur false-positive costs. The grader also accepts event tables with exactly `record_id,t`; the evaluator shows how to join them to the private scoring manifest. [C1 definition](QPD_LED_GENERALIZATION_SCORER_C1.txt)

For the **3 ms LED excess**, open `fixed-scores/pulse_candidate.json`. It pairs the independent validation arms at SNR 2, mean 1 versus 0, and uses every complete interval from a scheduled rise up to, but excluding, rise + 3 ms. Compare reconstructed and true active-minus-zero counts per pulse, bias, record-cluster standard error, raw excess SNR and FP per 100 pulses in each arm. Window FP uses the original whole-record C1 assignment. A claim at the native FP budget requires candidate FP/100 pulses ≤ native separately in both arms on identical exposure. Report both rates even when a budget fails. The display reuses replicate seeds across conditions, so this independent-arm calculation uses fixed validation.

## 4. Make the two sweep figures and the waveform inspection

```bash
python scripts/prepare_led_controlled.py --split display --replicates 64 --out "$D/display"
python scripts/evaluate_led_candidate.py --data "$D/display" --candidate "$A" \
    --label "My method v1" --out "$D/display-report"
python scripts/plot_led_nominal_diagnostics.py --data "$D/display-report" --out "$D/display-report/figures"
```

For an iteration comparison, add `--previous path/to/previous.py --previous-label "My method v0"` to the evaluation command. Keep each method's checkpoint and operating point fixed and identify them in the figure labels or report. Each one-second record goes through every method before any display cropping.

| Output in `figures/` (PNG and PDF) | What to inspect |
|---|---|
| `efficiency_nominal` | Strict one-to-one hits / all truth flips, versus LED mean at SNR 2 and versus SNR at mean 1 |
| `purity_nominal` | The same hits / all reconstructed flips. A point with zero predictions has undefined purity. |
| `waveform_nominal_0p5s` | Measured I and Q, LED rises, truth and one timestamp row per method. Five consecutive 100 ms panels show the entire first half-second. |

The sweeps pool counts across 64 records per point and show 95% record-bootstrap intervals. Strict matching maximizes the number of one-to-one pairs within ±200 µs at 10 kHz. C1 optimizes a different objective and can assign a late pair. Use the strict counts for these curves. The nominal display uses mean 1, SNR 2, identity I/Q augmentation, first onset 10 ms and replicate 0 (seed 106111000). Compare all methods on that same realization across iterations. The plot’s three-sample average is a visual aid.

Keep `plan.json`, `records.json`, `nominal_example.npz`, `figures/curves.json` and `figures/paired_nominal.json` with the C1 and pulse reports. The paired file also gives independently matched fixed-50-µs counts and transitions between unmatched, late and precise C1 truth states relative to native. Show the three plots inline when reporting a campaign iteration. For a newly trained SIDERIUS model, use the formal GPU training framework and include epoch-by-epoch training and fixed, record-disjoint validation loss, best and final loss/epochs, selected checkpoint and stopping reason. Label a frozen-checkpoint evaluation accordingly.
