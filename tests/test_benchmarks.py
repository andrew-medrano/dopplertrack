"""The benchmark configs load as events, and benchmarks/run.py scores a results table."""
import importlib.util
import json
from pathlib import Path

import pandas as pd
import pytest

from dopplertrack.tracks import EventFit

ROOT = Path(__file__).resolve().parents[1] / "benchmarks"
CONFIGS = sorted(ROOT.glob("*.json"))


def _run_module():
    spec = importlib.util.spec_from_file_location("bench_run", ROOT / "run.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_there_are_five_configs():
    assert {p.stem for p in CONFIGS} == {"oyster_bay", "eagle_cam", "railcam", "mcmaster", "safari_road"}


@pytest.mark.parametrize("path", CONFIGS, ids=lambda p: p.stem)
def test_config_events_load_and_coordinates_are_coarse_or_public(path):
    cfg = json.loads(path.read_text())
    evs = [EventFit(**e) for e in cfg["events"]] + [EventFit(**e) for e in cfg["extra"]]
    assert len(cfg["events"]) == 2 and len(cfg["harmonics"]) == 2
    assert all(e.t_start < e.t_end and e.f0 > 0 and e.d > 0 for e in evs)
    assert len(cfg["extra_harmonics"]) in (0, len(cfg["extra"]))
    ref = cfg["reference"]
    if path.stem not in ("oyster_bay", "eagle_cam"):       # the two public points upstream publishes
        assert round(ref["lat"], 2) == ref["lat"] and round(ref["lon"], 2) == ref["lon"]


def test_score_ranks_the_truth_pair_and_measures_distance(tmp_path):
    run = _run_module()
    cfg = json.loads((ROOT / "oyster_bay.json").read_text())
    ref = cfg["reference"]
    res = pd.DataFrame({
        "icao1": ["a1dfeb", "a1dfeb", "b00001", "b00002"], "icao2": ["ace391", "ace391", "b00003", "b00004"],
        "lat": [ref["lat"] + 0.001, ref["lat"], 30.0, 35.0], "lon": [ref["lon"], ref["lon"], -90.0, -100.0],
        "score": [30.0, 20.0, 27.0, 25.0]})
    res.to_csv(tmp_path / "results.csv", index=False)
    s = run._score(tmp_path / "results.csv", cfg)
    assert s["pairs"] == 3 and s["truth_rank"] == 1 and s["winner"] == "a1dfeb/ace391"
    assert abs(s["truth_km"] - 0.111) < 0.002 and s["hit"]
    assert abs(s["margin_pct"] - 100 * 3 / 27) < 0.1
