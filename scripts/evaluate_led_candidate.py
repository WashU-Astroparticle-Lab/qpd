"""Compare a frozen reconstruction adapter to native HMM-SF on identical records."""
import argparse
import hashlib
from importlib.util import module_from_spec, spec_from_file_location
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

from qpd.mlebench.led_generalization import write_json
from qpd.mlebench.led_generalization_hmmsf import predict_hmmsf
from qpd.mlebench.led_generalization_score import grade_generalization, score_record

CONFIG = Path(__file__).resolve().parents[1] / "configs/siderius/led_standalone_snr2.json"
TIMING = ("sample_rate_hz", "first_pulse_time_s", "on_time_s", "off_time_s")


def load_adapter(path):
    path = path.resolve()
    spec = spec_from_file_location("led_adapter_" + hashlib.sha256(str(path).encode()).hexdigest(), path)
    module = module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module.reconstruct


def event_table(rows, key):
    events = ((r["record_id"], t) for r in rows for t in
              (r["truth"] if key == "truth" else r["predictions"][key]))
    return pd.DataFrame(events, columns=["record_id", "t"])


def evaluate(data, out, candidate, label, previous, previous_label):
    receipt = json.loads((data / "COMPLETE.json").read_text())
    manifest = json.loads((data / "public/manifest.json").read_text())
    panel = json.loads((data / "private/panel.json").read_text())
    design = {r["record_id"]: r for r in json.loads((data / "private/design.json").read_text())}
    if len(manifest) != receipt["n_records"] or len(manifest) != len(design) or {r["record_id"] for r in manifest} != design.keys():
        raise ValueError("Incomplete panel or inconsistent record IDs")
    emg = json.loads(CONFIG.read_text())["baseline"]["nominal_emg"]
    adapters = {"candidate": load_adapter(candidate)}
    methods = [dict(key="hmmsf", label="Native HMM-SF", color="#27855d", marker="s"),
               dict(key="candidate", label=label, color="#c55622", marker="o")]
    if previous:
        adapters["candidate_previous"] = load_adapter(previous)
        methods.append(dict(key="candidate_previous", label=previous_label, color="#4775b9", marker="^"))
    out.mkdir(parents=True, exist_ok=False)
    rows, meta = [], []
    for row in manifest:
        rid = row["record_id"]
        timing = {k: row[k] for k in TIMING}
        with np.load(data / "public" / row["file"]) as f:
            iq = f["I"] + 1j * f["Q"]
        native = predict_hmmsf(iq, **timing, nominal_emg=emg).flip_times
        predictions = {"hmmsf": native.tolist()}
        for name, predict in adapters.items():
            values = np.asarray(predict(iq.copy(), **timing))
            if values.ndim != 1 or values.dtype.kind not in "iuf":
                raise ValueError(f"{name} must return a 1-D array of real timestamps for {rid}")
            predictions[name] = np.sort(values.astype(float)).tolist()
        with np.load(data / "private/truth" / f"{rid}.npz") as f:
            truth = f["flip_times"].tolist()
        d = design[rid]
        m, s = (d["parameters"][k] for k in ("mean_tunnels", "snr_max"))
        record = dict(record_id=rid, mean_tunnels=m, snr_max=s, truth=truth, predictions=predictions)
        meta.append(dict(**row, component=d["component"], zero_led=m == 0))
        if panel["split"] == "display":
            record["replicate"] = d["replicate"]
            record["metrics"] = {k: score_record(truth, v, timing["sample_rate_hz"], row["duration_s"])
                                 for k, v in predictions.items()}
            nominal = panel["nominal"]
            if (m, s, d["replicate"]) == (nominal["mean_tunnels"], nominal["snr_max"], 0):
                period = timing["on_time_s"] + timing["off_time_s"]
                first = timing["first_pulse_time_s"]
                count = int(np.ceil((row["duration_s"] - first) / period))
                pulses = first + period * np.arange(count)
                pulses = pulses[pulses < row["duration_s"]]
                np.savez_compressed(out / "nominal_example.npz", I=iq.real, Q=iq.imag,
                                    pulse_times=pulses, truth=truth, **predictions)
        rows.append(record)
    frame = pd.DataFrame(meta)
    tables = {k: event_table(rows, k) for k in ["truth", "hmmsf", *adapters]}
    for method in methods:
        key = method["key"]
        report = grade_generalization(tables[key], tables["truth"], frame, native_submission=tables["hmmsf"])
        write_json(out / f"c1_{key}.json", report)
    if panel["kind"] == "controlled" and panel["split"] == "validation":
        nominal = panel["nominal"]
        ids = {r["record_id"] for r in rows if r["snr_max"] == nominal["snr_max"]
               and r["mean_tunnels"] in (0, nominal["mean_tunnels"])}
        matched = frame[frame.record_id.isin(ids)].assign(pulse_condition="nominal")
        sub = {k: t[t.record_id.isin(ids)] for k, t in tables.items()}
        for key in adapters:
            result = grade_generalization(sub[key], sub["truth"], matched, native_submission=sub["hmmsf"])
            write_json(out / f"pulse_{key}.json", result["led_window_3ms"]["by_condition"]["nominal"])
    write_json(out / "plan.json", dict(**panel, plot_methods=methods,
               adapter_sha256={k: hashlib.sha256(p.read_bytes()).hexdigest() for k, p in
                               (("candidate", candidate), ("candidate_previous", previous)) if p},
               plot_note="Frozen operating points on identical full records."))
    write_json(out / "records.json", rows)
    write_json(out / "COMPLETE.json", dict(complete=True, records=len(rows)))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--data", required=True, type=Path)
    p.add_argument("--out", required=True, type=Path)
    p.add_argument("--candidate", required=True, type=Path, help="Python file defining reconstruct(iq, **timing)")
    p.add_argument("--label", default="Candidate", help="Method and checkpoint label for figures")
    p.add_argument("--previous", type=Path)
    p.add_argument("--previous-label", default="Previous")
    a = p.parse_args()
    evaluate(a.data, a.out, a.candidate, a.label, a.previous, a.previous_label)


if __name__ == "__main__":
    main()
