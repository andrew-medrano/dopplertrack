# Changelog

## 0.4.0 (2026-11-06)

First release of this fork ([andrew-medrano/dopplertrack](https://github.com/andrew-medrano/dopplertrack))
of [mjshiggins/dopplertrack](https://github.com/mjshiggins/dopplertrack).
It includes upstream's changes after 0.3.0 (listed at the end of this
entry). The search runs on a compact track store and scores in two stages,
with a 50 m refinement. On the Oyster Bay day the store is 646 MB against
6.2 GB for the 0.3 index, and a nationwide search took 88 s on a loaded M3 Pro
against 34 minutes for 0.3 on the same machine (swapping; 58 s against
309 s on the smaller eagle-cam day). docs/fast-search.md has the conditions.

### Breaking: searches run on a track store, in two stages

No compatibility layer: the 1 Hz day index, its loader and the
score-everything search are gone; `dopplertrack store` converts an old index
once.

| before | now |
|---|---|
| `fetch DATE day.npz` (1 Hz day index, 6.2 GB for a CONUS day) | `fetch DATE DIR` (track store, 646 MB for the same day) |
| — | `store DAY.npz DIR` converts a 0.3 day index (memory-mapped; 146 s and a 2.8 GB memory peak for the 6.2 GB Oyster Bay day, on a loaded machine) |
| `search … --day day.npz` scores every candidate | `search … --day DIR`: table scores for every candidate, exact scores for the best 20 points of the best 500 pairs, a 50 m refinement around the best 50 (`--exact-pairs`, `--exact-points`, `--refine-pairs`) |
| `refine … --day day.npz` | `refine … --day DIR` |
| — | `verdict RESULTS.csv` |

- `dopplertrack.store.TrackStore`: knots every 15 s with positions and
  velocities from a 31 s quadratic fit (cubic Hermite between knots), plus
  the 10 s light index the join and the silence filter use; float32 `.npy`
  files and a versioned `meta.json`, memory-mapped on load. It replaces
  `search.DayIndex`, `search.Segment` and `build_segments`. Exact scoring
  models curves from the stored tracks with their own velocities instead
  of differencing 1 Hz positions.
- `dopplertrack.cpa.closest_approach`: the closest approach to the stored
  curve, not its chords: the nearest chord among all chords of the window,
  then Newton steps on the arc of every chord within its arc bound. On
  turning and circling tracks it is within 0.3 m of brute force on the
  curve; a nearest-knot search had been off by up to 183 m there.
  `closest_approach_dense` is the brute-force reference.
- `score.EventWindow` scores on a log-frequency axis (0.1 % steps) with the
  harmonics summed once per window, in brightness and ridge mode alike.
  Scores correlate 0.999 with the linear-axis bilinear scorer on the
  synthetic tests; the best placement moves by at most one grid step.
  `locked_score`, `single_score` and `modelled_curve` take `store.Track`s.
- `dopplertrack.tables`: straight-pass score tables per event over distance ×
  speed × shift; `search.table_scores`, `shortlist` and `refine_candidates`
  are the stages around them. `run_search` reports counts and stage times in
  `.attrs["stats"]` and the CLI prints them. Results carry a `stage` column
  (`exact` or `refine`).
- The fine lattice's pre-gate judges distance only: on a turning track the
  chords' CPA time and speed can be off by more than any small slack, and
  a 3 s / 2 % slack dropped 2.3 % of the valid candidates of 20,000 pairs of
  the eagle-cam day.
- The silence filter takes candidate points in 10-minute groups, each against
  the flights near it in its own window, and checks only the points within
  reach of each flight: a pair whose flights stay together for an hour costs
  1 s of CPU instead of 22 s with one window for all its points.
- Terrain tiles are decoded once into memory-mapped float32 `.npy` files.
- Time shifts are swept every 0.25 s instead of 0.5 s (`score.SHIFTS`,
  `tables.DELTAS`). Event 2 takes event 1's best shift, so with 0.5 s steps a
  0.1 s change in a candidate's heard CPA time (the exact closest approach
  against the chords) could move event 1's pick by a step and event 2's
  score by half a point. On the railcam that dropped the true pair from 26th
  to 53rd after the exact stage, just outside the 50 refined pairs; with
  0.25 s steps it is 33rd there and 4th after refinement. The tables and the
  exact stage take about 2 s longer per search.
- `--utc` reads an ISO-8601 time without an offset as UTC, as its help says;
  it was read in the machine's time zone (an upstream bug), so on a machine
  in New York the Oyster Bay flyover's time matched no flight.
- The verdict leaves out every pair that shares a flight with the winner, in
  either role, as its docstring says; it compared flight 1 with flight 1 and
  flight 2 with flight 2 only. A one-event search (`--utc`) had an empty null
  and always came out "insufficient"; its null is now the other flights.
- The store's arrays are held as plain ndarray views of their memory maps:
  slicing an `np.memmap` builds a new memmap object each time, which cost
  6 % of the lens stage.
- On the synthetic scene of tests/test_search.py the two-stage search puts
  the winner 2 m from the observer; scoring every 200 m candidate put it
  143 m away. On real data (0.3 indexes converted with `store`): Oyster Bay
  ranks the true pair 1 of 500 shortlisted pairs, 26 m from the answer, 8.1 %
  ahead, located (p = 0.09; 0.05 with its extra flyover); the eagle cam
  ranks it 13 of 500 with two events and 1 of 200 with its three extra
  flyovers (4.9 % ahead, ambiguous, p = 0.42). The railcam's true pair
  ranks 53 of 499 (the prototype had it 3rd). McMaster Forest, held out, is
  located 0.76 km from its reference point with the silence filter off and
  its extras (p = 0.001), as in the prototype; Safari Road is missed.
- The silence filter skips points it has already removed: on the McMaster
  day with the filter on, 1,181 s of CPU instead of 2,046 s, same result.
- Against 0.3 on the same machine, the eagle-cam search is 5.4× faster in
  wall time (309 s to 58 s) and 5.8× in CPU time; Oyster Bay 23× and 18×,
  inflated by 0.3 swapping. The verdict is weaker: with the extras the
  eagle cam is ambiguous (p = 0.42) where 0.3 said located (p = 0.05),
  because the null is the shortlist.

### Documentation and repository

- README: a fork header with credits, install and a worked example for the
  new commands, measured performance, how the fast search works, responsible
  use, and data licences.
- docs/fast-search.md: the design, storage and timing of the fast search,
  with the prototype measurements and what did not work.
- docs/heldout.md: the held-out protocol and the two held-out recordings
  (McMaster Forest: one hit, with the silence filter off after the
  pre-registered run missed; Safari Road: missed in every variant).
- docs/other-signals.md: what else the recordings carry (jet loudness, mains
  hum, rest frequency, stereo), as leads.
- benchmarks/: JSON configs for the five recordings (events, harmonics, truth
  pairs, UTC days, audio sources; site coordinates rounded to 0.01°) and
  `benchmarks/run.py`, which runs a search and scores it.
- pyproject: version 0.4.0, Andrew Medrano added as an author, URLs point to
  the fork. The release workflow no longer publishes to PyPI (the job is
  kept, disabled), since the fork would publish under upstream's name.
- .gitignore: audio, arrays, stores, parquet and pickle files and the research
  work directories.

### From upstream after 0.3.0

- `dopplertrack.verdict`: the winner of a search is judged against the wrong
  pairs the same search scored (an exponential upper tail fitted on the best
  wrong pairs; p is the chance a wrong winner would lead the next wrong pair by
  the observed gap). `search` prints it after the top table: located, ambiguous,
  no match, or insufficient. On the saved results two events alone are always
  ambiguous, the extra flyovers make both recordings located, and the box moved
  150 km from the answer is no match. Those tables are now under docs/results
  and tests/test_verdict.py checks them.

- `detect` flags events that coincide with a broadband burst (`burst`: the
  fraction of the 3 min around the CPA in which the whole band is loud at
  once, three robust sigmas above its median level). On a railcam recording
  every top detection was a locomotive passing the camera, whose engine tone
  sweeps like a low-altitude flyover and arrives with the rumble.
  `--max-burst` drops them.

- A third ground truth, a 12-hour Railstream camera in Atlanta, searched
  blind: negative. The true pairs rank 6–32 nationwide within 440 m and never
  win; the archived stream has a 107 s dropout that puts events on two clocks;
  non-aircraft extras bury the right pair; and two of twenty verdicts with no
  right answer were "located" (docs/experiments.md). Its tables are in
  docs/results and the two false positives are expected failures in
  tests/test_verdict.py.

## 0.3.0 (2026-11-04)

The search no longer needs a box. One day of flights over a continent is
indexed once; flight pairs are joined on the geometry of the two events, ground
points are laid only in each pair's lens, and the same exact CPA gates, silence
filter and locked scorer rank them. Validated on both ground truths over the
whole continental US: Oyster Bay ranks the true pair 1 of 9,669 (32 m); the
eagle cam needs its further detected flyovers (`--extra`) and then ranks 1 of
13,298 (627 m). See docs/experiments.md.

**Breaking: the CLI is six commands.**

| before | now |
|---|---|
| `fetch-day DATE traces.npz --bbox …`, `extract-region DIR traces.npz --bbox …`, `index-day traces.npz day.npz` | `fetch DATE day.npz [--bbox …] [--archive DIR] [--traces traces.npz]`, which writes the day index directly |
| `search … --traces traces.npz --bbox …` (200 m grid over a box) | `search … --day day.npz [--bbox …]` (the join; `--bbox` restricts it) |
| `search … --point LAT LON` | `search … --day day.npz --point LAT LON` |
| `search-wide … --index day.npz` | `search … --day day.npz` |
| `search-single … --event EV --utc T` | `search … --events EV --utc T` |
| `refine … --traces traces.npz --bbox …` | `refine … --day day.npz` |

**Breaking: the library.** `dopplertrack.wide` is merged into
`dopplertrack.search` (`DayIndex`, `pair_join`, `single_join`, `lens_candidates`,
`score_candidates`, `run_search`, `rescore_with_extra`, `pick_flight`). The grid
search (`cpa_table`, `pair_candidates`, the CPA cache, `Grid.points`) and
`pipeline.run_search`, `run_single_search`, `extract_region`, `fetch_region`,
`segments_for_search`, `relocate_row` are gone; `pipeline.fetch` builds a day
index, `pipeline.refine` takes one. `build_segments` lost `clip_enu`. Result
CSVs carry `seg1`/`seg2` into the index they were made from, and `refine` checks
them against it.

Also:

- The locked scorer samples the spectrogram window along time once per event
  and walks the sweep in small chunks: twice as fast on one core, eleven times
  on sixteen (it was memory-bandwidth-bound in a process pool).
- The lens stage computes exact CPA with a KD-tree nearest-sample query instead
  of a dense distance matrix (12 min → 71 s for a CONUS day).
- `load_traces` reads the metadata member once (371 s → 5 s on 35k traces).
- Terrain tiles are written atomically and re-fetched if corrupt, so a process
  pool can share the cache.
- One BLAS thread per process is set before numpy loads in the CLI.

## 0.2.0 (2026-11-03)

- Coverage-weighted scoring: brightness summed over the whole measured window,
  so a flight that leaves the data mid-pass is not judged on its bright half.
- `refine --extra`: further detected flyovers scored on the same clock.
- `detect`: automatic Doppler event detection over a whole recording.
- `fetch-day`: stream an adsb.lol day archive from GitHub without storing it.
- Eagle-cam validation (second ground truth).

## 0.1.0 (2026-11-03)

- Spectrogram, curve fit, regional grid search, refine; Oyster Bay result.
