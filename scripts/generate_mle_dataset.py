#!/usr/bin/env python
"""CLI to generate the MLE-bench QPD datasets (signal + background).

Examples
--------
Plan only (no data written)::

    python scripts/generate_mle_dataset.py --dry-run

Full run, 2 GB / 1 h budget (whichever binds first)::

    python scripts/generate_mle_dataset.py --out-dir data/mle_bench

Tiny smoke test (2 short chunks/dataset into a temp dir)::

    python scripts/generate_mle_dataset.py --smoke
"""

from __future__ import annotations

import argparse
import json
import os
import sys

from qpd.mlebench.config import DEFAULT_PHYSICS, GenerationBudget
from qpd.mlebench.generate import generate_all, plan_generation


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--out-dir", default="data/mle_bench",
                   help="output root (one subdir per competition)")
    p.add_argument("--max-bytes", type=float, default=2e9,
                   help="total byte budget across both datasets (default 2e9)")
    p.add_argument("--max-compute-seconds", type=float, default=3600.0,
                   help="total wall-clock budget across both datasets (s)")
    p.add_argument("--max-chunks-per-dataset", type=int, default=None,
                   help="hard cap on chunks per dataset (safety net)")
    p.add_argument("--chunk-seconds", type=float, default=30.0)
    p.add_argument("--train-fraction", type=float, default=0.6)
    p.add_argument("--seed", type=int, default=0, help="master seed")
    p.add_argument("--sec-per-chunk", type=float, default=1.0,
                   help="compute-cost estimate used for budget planning")
    p.add_argument("--dry-run", action="store_true",
                   help="print the generation plan and exit")
    p.add_argument("--smoke", action="store_true",
                   help="tiny end-to-end test: 2 chunks/dataset, 2 s each")
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)

    if args.smoke:
        import tempfile
        from qpd.mlebench.generate import generate_all as _gen
        budget = GenerationBudget(
            max_total_bytes=1e9, max_total_compute_seconds=600.0,
            max_chunks_per_dataset=2,
        )
        out = args.out_dir if args.out_dir != "data/mle_bench" else (
            tempfile.mkdtemp(prefix="qpd_mle_smoke_"))
        print(f"[smoke] writing to {out}")
        res = _gen(out, budget=budget, chunk_seconds=2.0,
                   train_fraction=0.5, master_seed=args.seed)
        print(json.dumps(res, indent=2, default=str))
        return 0

    budget = GenerationBudget(
        max_total_bytes=args.max_bytes,
        max_total_compute_seconds=args.max_compute_seconds,
        max_chunks_per_dataset=args.max_chunks_per_dataset,
    )

    plan = plan_generation(
        budget,
        chunk_seconds=args.chunk_seconds,
        sample_rate_hz=DEFAULT_PHYSICS.sample_rate_hz,
        n_datasets=2,
        sec_per_chunk=args.sec_per_chunk,
    )
    print("Generation plan:")
    print(json.dumps(plan, indent=2))

    if args.dry_run:
        return 0

    os.makedirs(args.out_dir, exist_ok=True)
    result = generate_all(
        args.out_dir,
        budget=budget,
        chunk_seconds=args.chunk_seconds,
        train_fraction=args.train_fraction,
        master_seed=args.seed,
        sec_per_chunk=args.sec_per_chunk,
    )
    print("\nDone. Summary:")
    print(json.dumps(result, indent=2, default=str))
    return 0


if __name__ == "__main__":
    sys.exit(main())
