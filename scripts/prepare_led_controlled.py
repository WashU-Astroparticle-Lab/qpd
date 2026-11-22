"""Generate fixed-condition training, validation, or recurring display waveforms."""
import argparse
from pathlib import Path
from uuid import UUID

import numpy as np

from qpd.mlebench.led_dataset import LEDDatasetConfig, build_led_context
from qpd.mlebench.led_generalization import load_spec, simulate_record, write_json

CONFIG = Path(__file__).resolve().parents[1] / "configs/siderius/led_standalone_snr2.json"
MEANS = [0., .1, .3, .5, 1., 2., 3., 5., 10.]
SNRS = [.5, .75, 1., 1.5, 2., 3., 4., 5.]
SPLIT = {"train": 1, "validation": 2, "display": 3}
DISPLAY_SEED = 106111000


def generate(split, replicates, cross, out):
    spec = load_spec(CONFIG)
    nominal = {key: v["nominal"] for key, v in spec["parameters"].items()}
    nominal.update(first_pulse_time_s=.010, burst_mu_s=spec["burst_mu_s"])
    mean, snr = nominal["mean_tunnels"], nominal["snr_max"]
    conditions = ([(m, s) for m in MEANS for s in SNRS] if cross else
                  [(m, snr) for m in MEANS] + [(mean, s) for s in SNRS if s != snr])
    out.mkdir(parents=True, exist_ok=False)
    for folder in ("public/waveforms", "private/truth"):
        (out / folder).mkdir(parents=True)
    cfg = LEDDatasetConfig(**{k: v for k, v in nominal.items() if k not in ("mean_tunnels", "snr_max")})
    context = build_led_context(cfg)
    public, design = [], []
    for m, s in conditions:
        params = dict(nominal, mean_tunnels=m, snr_max=s)
        for rep in range(replicates):
            key = [114902, SPLIT[split], MEANS.index(m), SNRS.index(s), rep]
            id_rng = np.random.default_rng(np.random.SeedSequence(key + [0]))
            rid = UUID(bytes=id_rng.bytes(16)).hex
            seed = (DISPLAY_SEED + rep if split == "display" else
                    int(np.random.default_rng(np.random.SeedSequence(key + [1])).integers(2**63-1)))
            result = simulate_record(context, params, cfg.duration_s, seed)
            iq = result.iq.astype(np.complex64)
            file = f"waveforms/{rid}.npz"
            np.savez_compressed(out / "public" / file, I=iq.real, Q=iq.imag)
            np.savez_compressed(out / "private/truth" / f"{rid}.npz", flip_times=result.flip_times)
            public.append(dict(record_id=rid, file=file, n_samples=len(iq), duration_s=cfg.duration_s,
                               sample_rate_hz=cfg.sample_rate_hz, first_pulse_time_s=cfg.first_pulse_time_s,
                               on_time_s=cfg.on_time_s, off_time_s=cfg.off_time_s))
            design.append(dict(record_id=rid, component="nominal", parameters=params,
                               simulation_seed=seed, replicate=rep))
    np.random.default_rng(90432).shuffle(public)
    write_json(out / "public/manifest.json", public)
    write_json(out / "private/design.json", design)
    write_json(out / "private/panel.json", dict(kind="controlled", split=split, replicates=replicates,
               nominal=nominal, mean_grid=MEANS, snr_grid=SNRS, conditions=conditions,
               example_replicate=0, seed_base=DISPLAY_SEED if split == "display" else None))
    write_json(out / "COMPLETE.json", dict(n_records=len(public)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--split", choices=SPLIT, required=True)
    parser.add_argument("--replicates", type=int, default=64, help="Independent records per setting within a split")
    parser.add_argument("--cross", action="store_true", help="Include all 72 mean/SNR combinations")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if args.replicates < 1 or (args.cross and args.split == "display"):
        parser.error("Choose a positive replicate count and use --cross only with train or validation")
    generate(args.split, args.replicates, args.cross, args.out)


if __name__ == "__main__":
    main()
