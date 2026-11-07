"""The command line end to end on the synthetic scene of test_search: convert a 0.3 day index
with `store`, then `search` (two events, a known point, one event at a known time), `refine`
and `verdict` on the files the previous step wrote."""
import json

import pandas as pd
import pytest

from dopplertrack import cli, tracks
from dopplertrack.geo import haversine_m
from test_search import C, LAT0, LON0, T0, _audio, _flights, _index
from test_store import _day_index_npz


@pytest.fixture(scope="module")
def files(tmp_path_factory):
    d = tmp_path_factory.mktemp("cli")
    st = _index()
    sp = _audio(st, (("aaa001", 120.0), ("aaa002", 130.0)))
    sp.save(d / "rec.npz")
    for name, win, band in (("ev1", (230, 400), (90, 170)), ("ev2", (1150, 1310), (100, 180))):
        ev = tracks.fit_event(sp, *win, *band, c=C, objective="brightness")
        (d / f"{name}.json").write_text(json.dumps(ev.to_dict()))
    _day_index_npz(d / "day.npz", _flights())
    return d


def test_store_then_search_refine_verdict(files, capsys):
    d = files
    cli.main(["store", str(d / "day.npz"), str(d / "store")])
    assert (d / "store" / "meta.json").exists()
    cli.main(["search", str(d / "rec.npz"), "--events", str(d / "ev1.json"), str(d / "ev2.json"), "--day", str(d / "store"),
              "--out", str(d / "res.csv"), "--no-dem", "--workers", "1"])
    out = capsys.readouterr().out
    assert "stage times" in out and "verdict:" in out
    res = pd.read_csv(d / "res.csv", keep_default_na=False)
    assert set(res.stage) == {"exact", "refine"}
    best = res.iloc[0]
    assert (best.icao1, best.icao2) == ("aaa001", "aaa002") and haversine_m(best.lat, best.lon, LAT0, LON0) < 150.0
    cli.main(["refine", str(d / "rec.npz"), str(d / "res.csv"), "--events", str(d / "ev1.json"), str(d / "ev2.json"),
              "--day", str(d / "store"), "--out", str(d / "fine.csv"), "--no-dem", "--half", "200", "--spacing", "100"])
    fine = pd.read_csv(d / "fine.csv")
    assert len(fine) >= 16 and haversine_m(fine.lat[0], fine.lon[0], LAT0, LON0) < 250.0
    capsys.readouterr()
    cli.main(["verdict", str(d / "res.csv")])
    assert capsys.readouterr().out.startswith("verdict: insufficient")       # two pairs is no crowd


def test_search_known_point_and_known_time(files):
    d = files
    if not (d / "store" / "meta.json").exists():
        cli.main(["store", str(d / "day.npz"), str(d / "store")])
    cli.main(["search", str(d / "rec.npz"), "--events", str(d / "ev1.json"), str(d / "ev2.json"), "--day", str(d / "store"),
              "--out", str(d / "pt.csv"), "--no-dem", "--workers", "1", "--point", str(LAT0), str(LON0)])
    pt = pd.read_csv(d / "pt.csv", keep_default_na=False)
    assert len(pt) == 1 and (pt.icao1[0], pt.icao2[0]) == ("aaa001", "aaa002") and set(pt.stage) == {"exact"}
    ev1 = json.loads((d / "ev1.json").read_text())
    cli.main(["search", str(d / "rec.npz"), "--events", str(d / "ev1.json"), "--utc", str(T0 + ev1["t_cpa"]),
              "--tol", "5", "--day", str(d / "store"), "--out", str(d / "one.csv"), "--no-dem", "--workers", "1"])
    one = pd.read_csv(d / "one.csv", keep_default_na=False)
    assert len(one) and (one.seg2 == -1).all() and one.icao1[0] in ("aaa001", "ccc001")


def test_search_rejects_a_day_index(files):
    with pytest.raises(ValueError, match="dopplertrack store"):
        cli.main(["search", str(files / "rec.npz"), "--events", str(files / "ev1.json"), str(files / "ev2.json"),
                  "--day", str(files / "day.npz"), "--out", str(files / "x.csv"), "--no-dem"])
