"""Plan or generate the campaign's mixed train and validation splits."""
import argparse
from pathlib import Path

import numpy as np

from qpd.mlebench.led_generalization import (
    generate_pilot, load_spec, plan_pilot, plan_summary, write_json,
)

ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "configs/siderius/led_standalone_snr2.json"
SIZES = {"train": 20000, "validation": 4000}
SEEDS = {"train": 410010, "validation": 510010}


def in_holdout(p):
    return ((p["sample_rate_hz"] <= 10000 and p["snr_max"] <= 1.4 and p["burst_tau_s"] >= .00085)
            or (p["sample_rate_hz"] >= 100000 and p["off_time_s"] <= .010 and p["burst_tau_s"] >= .001))


def plan_scientific(spec, split, count):
    """Replay the scientific design in shuffled blocks with exact six-stratum quotas."""
    if count < 200 or count % 200 or count > SIZES[split]:
        raise ValueError(f"Use a multiple of 200 up to {SIZES[split]} for {split}")
    base = SEEDS[split]
    blocks = count // 200
    unit = {(c, zero): n * (1 if zero else 9) // 10
            for c, n in [("broad", 120), ("nominal", 50), ("stress", 30)]
            for zero in (False, True)}
    buckets = {key: [] for key in unit}
    for batch in range(10000):
        for row in plan_pilot(spec, n_records=2000, seed=base + batch):
            p = row["parameters"]
            key = row["component"], p["mean_tunnels"] == 0
            if len(buckets[key]) < blocks * unit[key] and not in_holdout(p):
                row.update(scientific_split=split, evaluation_partition="iid", heldout_cell=None)
                buckets[key].append(row)
        if all(len(rows) == blocks * unit[key] for key, rows in buckets.items()):
            break
    else:
        raise RuntimeError("Unable to fill scientific design quotas")
    result = []
    for i in range(blocks):
        block = [row for key, n in unit.items() for row in buckets[key][i*n:(i+1)*n]]
        np.random.default_rng(np.random.SeedSequence([base, 0, i, 99])).shuffle(block)
        result.extend(block)
    for field in ("record_id", "simulation_seed"):
        if len({row[field] for row in result}) != count:
            raise ValueError(f"Duplicate {field}")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--split", required=True, choices=SIZES)
    parser.add_argument("--count", type=int, help="Stratified prefix, in multiples of 200")
    parser.add_argument("--out", type=Path, help="Generate waveforms in a new directory")
    args = parser.parse_args()
    spec = load_spec(CONFIG)
    rows = plan_scientific(spec, args.split, SIZES[args.split] if args.count is None else args.count)
    print(plan_summary(rows), flush=True)
    if args.out:
        sources = [CONFIG, Path(__file__)] + sorted((ROOT / "src/qpd").rglob("*.py"))
        generate_pilot(spec, rows, args.out, source_paths=sources, scientific_split=args.split)
        write_json(args.out / "private/panel.json", {"kind": "mixed", "split": args.split})


if __name__ == "__main__":
    main()
