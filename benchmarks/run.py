"""Run one benchmark search and score it against its reference.

    uv run python benchmarks/run.py benchmarks/oyster_bay.json --spec rec.npz --day store/ --out runs/oyster
    uv run python benchmarks/run.py benchmarks/oyster_bay.json --score-only --out runs/oyster

Writes the config's events as event JSON into OUT, runs `dopplertrack search` with the
config's harmonics, extras and search_args (plus any arguments after `--`), then prints, for
the two-event results and the results with extras: the winner, its distance from the
reference point, the rank of the truth pair (when the config names one) and the verdict.
A summary goes to OUT/score.json.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

import pandas as pd

from dopplertrack.geo import haversine_m
from dopplertrack.verdict import best_per_pair, verdict


def _score(csv: Path, cfg: dict) -> dict:
    res = pd.read_csv(csv, dtype={"icao1": str, "icao2": str})
    b = best_per_pair(res).reset_index(drop=True)
    ref = cfg["reference"]
    d = haversine_m(b.lat.to_numpy(), b.lon.to_numpy(), ref["lat"], ref["lon"])
    w = b.iloc[0]
    v = verdict(res)
    out = dict(pairs=len(b), winner=f"{w.icao1}/{w.icao2}", winner_km=round(float(d[0]) / 1000, 2),
               winner_score=round(float(w.score), 2), verdict=v.verdict, p=round(float(v.p), 3),
               hit=bool(d[0] <= 2000))
    if len(b) > 1:
        out["margin_pct"] = round(100 * float((b.score.iloc[0] - b.score.iloc[1]) / b.score.iloc[1]), 1)
    tp = cfg.get("truth_pair")
    if tp:
        m = ((b.icao1 == tp[0]) & (b.icao2 == tp[1])).to_numpy().nonzero()[0]
        out["truth_rank"] = int(m[0]) + 1 if len(m) else None
        out["truth_km"] = round(float(d[m[0]]) / 1000, 3) if len(m) else None
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("config")
    ap.add_argument("--spec", help="spectrogram .npz of the benchmark's audio")
    ap.add_argument("--day", help="track store of the benchmark's UTC day")
    ap.add_argument("--out", required=True, help="directory for events, results and score.json")
    ap.add_argument("--score-only", action="store_true", help="only score existing results in OUT")
    a, rest = ap.parse_known_args(argv)
    cfg = json.loads(Path(a.config).read_text())
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    res, res_x = out / "results.csv", out / "results_extra.csv"
    if not a.score_only:
        if not (a.spec and a.day):
            ap.error("--spec and --day are needed unless --score-only")
        evs = []
        for i, e in enumerate(cfg["events"], 1):
            p = out / f"ev{i}.json"; p.write_text(json.dumps(e)); evs.append(str(p))
        xs = []
        for i, e in enumerate(cfg.get("extra", []), 1):
            p = out / f"x{i}.json"; p.write_text(json.dumps(e)); xs.append(str(p))
        cmd = [sys.executable, "-m", "dopplertrack.cli", "search", a.spec, "--events", *evs, "--day", a.day,
               "--harmonics", *map(str, cfg["harmonics"]), "--out", str(res), *cfg.get("search_args", [])]
        if xs:
            cmd += ["--extra", *xs, "--out-extra", str(res_x)]
            if cfg.get("extra_harmonics"):
                cmd += ["--extra-harmonics", *map(str, cfg["extra_harmonics"])]
        rest = [r for r in rest if r != "--"]
        print(" ".join(cmd + rest), flush=True)
        subprocess.run(cmd + rest, check=True)
    summary = {"two events": _score(res, cfg)}
    if res_x.exists():
        summary["with extras"] = _score(res_x, cfg)
    (out / "score.json").write_text(json.dumps(summary, indent=1))
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
