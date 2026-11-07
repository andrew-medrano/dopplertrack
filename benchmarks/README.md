# Benchmarks

Five recordings with a known or pre-registered answer, as small JSON configs,
and a script that runs the search on one and scores it. No audio, ADS-B or
track stores are included: the audio is linked, and the ADS-B day is
fetched from adsb.lol.

| config | recording | UTC day | role |
|---|---|---|---|
| `oyster_bay.json` | colsto's one-hour reference clip, Oyster Bay, NY | 2026-09-08 | development |
| `eagle_cam.json` | Southwest Florida Eagle Cam, 12 h | 2024-12-28 | development |
| `railcam.json` | Railstream Howell Wye camera, Atlanta, 11.9 h | 2026-05-28 | development |
| `mcmaster.json` | Hamilton Naturalists' Club, McMaster Forest, 12 h | 2026-05-27 | held out |
| `safari_road.json` | Hamilton Naturalists' Club, Safari Road, 7 h | 2023-07-09 | held out |

Each config holds the two events and any extra flyovers as `fit`/`detect`
output (times in seconds from the first sample of the audio as linked), the
harmonics per event, the extra search arguments, the truth pair's ICAO
addresses where one is known, and a reference point. Reference points of
the held-out sites and the railcam are rounded to 0.01° (about 1 km); the
held-out recorder positions are not known more precisely than that anyway.
How each held-out recording was assembled from its files is in its `audio`
note and in [docs/heldout.md](../docs/heldout.md).

## Running one

```bash
uv sync
# the audio, once (yt-dlp or a browser download; the configs name the source)
uv run dopplertrack spectrogram oyster_bay.wav runs/oyster_bay/rec.npz
# the ADS-B day, once: about 4 GB streamed, a track store of a few hundred MB written
uv run dopplertrack fetch 2026-09-08 stores/2026-09-08 --bbox 24 50 -125 -66
# search and score; anything after the known options goes to `dopplertrack search`
OBJC_DISABLE_INITIALIZE_FORK_SAFETY=YES uv run python benchmarks/run.py benchmarks/oyster_bay.json \
    --spec runs/oyster_bay/rec.npz --day stores/2026-09-08 --out runs/oyster_bay --workers 11
```

`run.py` writes the events to `OUT/ev*.json` and `OUT/x*.json`, runs
`dopplertrack search` with them, and prints the winner, its distance from
the reference point, the rank and distance of the truth pair, the margin
over the runner-up and the verdict, for the two events alone and with the
extras. The same summary goes to `OUT/score.json`. `--score-only` rescores
results already in `OUT`. Set `OBJC_DISABLE_INITIALIZE_FORK_SAFETY=YES` on
macOS only.

The McMaster and Safari Road configs carry `--no-silence`, the setting of
the second held-out run. The pre-registered first run had the silence
filter on; to repeat it, edit `search_args` to `[]`. The railcam config also
carries `--no-silence`, as upstream's railcam experiments did.

A day fetched with dopplertrack 0.3 (`DAY.npz`) converts with
`dopplertrack store DAY.npz DIR` instead of fetching again.

## Results

Scores from `run.py --score-only` for searches run with these configs on
2026-11-06 (this package, M3 Pro, 11 workers,
a loaded machine; times and conditions in
[docs/fast-search.md](../docs/fast-search.md)). The 0.3 day indexes of the
research were converted with `dopplertrack store` instead of fetched again.

| config | two events | with extras |
|---|---|---|
| oyster_bay | truth 1st, 26 m; located, p = 0.09 | truth 1st, 96 m; located, p = 0.05 |
| eagle_cam | truth 13th, 0.80 km; winner 203 km off; no match | truth 1st, 0.74 km; ambiguous, p = 0.42 |
| railcam | truth 53rd, 0.65 km; winner 691 km off; ambiguous | — |
| mcmaster (`--no-silence`) | winner 542 km off; no match | winner 1.19 km from the rounded point (0.76 km from the unrounded one); located, p = 0.001 |
| mcmaster, silence on | winner 284 km off; ambiguous | winner 3,492 km off; ambiguous |
| safari_road (`--no-silence`) | winner 1,715 km off; no match | winner 3,290 km off; no match |
| safari_road, silence on | winner 2,629 km off; no match | winner 3,450 km off; no match |

Distances to a rounded reference point can be off by up to about 0.7 km,
which matters near the 2 km hit line; the distances in docs/heldout.md use
the unrounded points.
