# Held-out recordings

The development recordings (Oyster Bay, the eagle cam, the railcam) were used
to build and tune the search, so their results are not a fair test. This page
covers two recordings whose answers were written down before any of their
audio was analysed, and the protocol they were run under.

Short version: two held-out sites, one hit. The pre-registered run missed
both. The hit (McMaster Forest, 0.7 km) came from a second run with the
silence filter off, a setting chosen after the first run missed. The other
site (Safari Road) missed in every variant, and its verdicts said so.

## The recordings

Both are long soundscape recordings published by the Hamilton Naturalists'
Club on the Internet Archive. The item metadata names the club as creator
and gives no licence, so the audio is not redistributed here.

| | McMaster Forest | Safari Road |
|---|---|---|
| Internet Archive item | [`HNCSW720260527`](https://archive.org/details/HNCSW720260527) | [`HNCAM01820230708`](https://archive.org/details/HNCAM01820230708) |
| place | McMaster Forest Nature Preserve, Ancaster, Ontario | a trail site on Safari Road, Flamborough, Ontario |
| reference point (rounded) | 43.24, -79.95 | 43.35, -80.12 |
| recorder position | unknown within the preserve | unknown within the site |
| window | 2026-05-27, 12:00-24:00 UTC (12 h, hourly files) | 2023-07-09, 10:00-17:02 UTC (7 h, 30 min files) |

Each truth file recorded the reference point with an 800 m radius, the
inputs allowed to the search (the audio and the UTC date only), and a
timestamp before the first download.

## Protocol v1 (pre-registered)

Written before any held-out audio was analysed. The truth was read only by
the final scoring step.

1. Place each file at its filename start time, silence in the gaps, decoded
   to 4 kHz mono.
2. `dopplertrack spectrogram` with defaults.
3. `dopplertrack detect --band 60 340 --top 15 --max-burst 0.15` over the
   whole window.
4. Events: the two highest-ranked detections, earlier one first. Extras: the
   next three detections.
5. Search the CONUS day of that UTC date, two events and then with the
   extras, in two ways: 0.3's search (log-axis scorer, every candidate scored
   exactly) and the research prototype of the fast search (500 pairs, 20
   points, 50 m refinement).
6. Score: distance from the winner to the reference point; a **hit** if
   within 2 km, a **near miss** if within 20 km.

Any rerun with changed settings is reported separately as such.

## Results

Distance from the winner to the reference point, and the verdict:

| run | McMaster Forest | Safari Road |
|---|---|---|
| v1, full search, 2 events | 864 km, no match (p = 0.98) | 3,388 km, no match (p = 1.00) |
| v1, full search, with extras | 886 km, no match (p = 0.82) | 1,696 km, no match (p = 0.75) |
| v1, fast, 2 events | 285 km, ambiguous (p = 0.49) | 2,012 km, no match (p = 0.53) |
| v1, fast, with extras | 3,492 km, ambiguous (p = 0.40) | 3,449 km, ambiguous (p = 0.40) |
| silence off, full search, 2 events | 542 km, no match (p = 0.90); the best pair within 800 m ranks 2nd | 1,918 km, ambiguous (p = 0.26) |
| silence off, full search, with extras | **0.68 km, located (p = 0.005)** | 743 km, ambiguous (p = 0.26) |
| silence off, fast, 2 events | 543 km, no match (p = 0.51); the best pair within 800 m ranks 2nd | 1,715 km, no match (p = 0.67) |
| silence off, fast, with extras | **0.76 km, located (p = 0.001)** | 1,715 km, no match (p = 0.98) |

The table above is from 2026-11-05: 0.3's search and the research prototype
of the fast search. The package in this repository, with the same events
(`benchmarks/mcmaster.json`, `benchmarks/safari_road.json`). The silence-off
rows are from the current code (2026-11-07: 0.25 s time-shift steps, and a
verdict that leaves out every pair sharing a flight with the winner); the
silence-on rows were not rerun and are from 2026-11-06, with 0.5 s steps:

| run (package) | McMaster Forest | Safari Road |
|---|---|---|
| silence on, 2 events | 284 km, ambiguous (p = 0.49) | 2,629 km, no match (p = 0.65) |
| silence on, with extras | 3,492 km, ambiguous (p = 0.42) | 3,450 km, no match (p = 0.86) |
| silence off, 2 events | 543 km, ambiguous (p = 0.49); the best pair within 800 m ranks 2nd (0.74 km) | 1,715 km, no match (p = 0.56) |
| silence off, with extras | **0.71 km, located (p = 0.0005)**, 23 % ahead | 1,715 km, no match (p = 0.84) |

With 0.5 s steps the silence-off McMaster results were 542 km, no match
(p = 0.53) for two events and 0.76 km, located (p = 0.001), 21 % ahead with
the extras; Safari Road's with extras was 3,290 km, no match (p = 0.63).

The package reproduces the prototype's held-out results: the same hit at
McMaster with the silence filter off and the extras, and misses otherwise.


**McMaster Forest.** The v1 run missed everywhere. Its two events are
5.3 hours apart, and the silence filter rejects a point if any other
propeller flight passed within 6 km between them. Upstream's railcam
experiments already turned the filter off on a 12-hour recording for this
reason ([experiments.md](experiments.md)). So a second run with
`--no-silence` was declared, after the v1 miss was seen but before the
second run's result was. With two events
alone it still missed, with the best pair within 800 m ranked 2nd. With the three extra
flyovers it located the recording 0.68 km from the reference point (inside
the 800 m radius). The winning pair is a Beechcraft Bonanza and a Pilatus
PC-12. The extras matched two more aircraft and a second pass of the same
Bonanza 4.6 hours later.

This is one hit, with a setting chosen after a miss. It is not a blind
result, and the reference point is the preserve, not the recorder.

**Safari Road.** Missed in every variant. Its detections are weak: the
strongest scores 8.0, against 17.2 at McMaster. The verdicts never said
"located" for a wrong answer here.

**The recorder clock.** At McMaster the winning pair puts the audio 593 s
later than the filenames say: the recorder's clock was about 10 minutes
behind UTC. The search only uses time differences between events, so this
did not matter. Never assume a recording's clock.

## Protocol v2: every detection, no silence filter

Written after McMaster's v1 result and its diagnostic, and before Safari
Road's v1 score was read. It uses the top 6 detections, tries all 15 pairs
of them as seeds, and scores each candidate location and clock against every
detection, with each one matched to its best flight or to nothing.

It was tried first on the Oyster Bay development recording, where the
answer is known. It missed: the best candidate by exact score was 136 km
off, and by table score 1,552 km off, after 609 s. A voting variant also
missed. Since it failed where 0.3's search works, it was not run on the
held-out sites.

## Fit quality

On Oyster Bay, `detect`'s events were compared with the real tracks of the
flights they matched. On clean straight passes the automatic fit is within
1-3 % of the truth in speed and distance (the PC-12: speed +1 %, distance
+3 %, CPA time -0.4 s). On manoeuvring or faint passes it is badly wrong:
speed -56 % to +63 %, distance down to -94 %. Six of `detect`'s top ten on
Oyster Bay are not aircraft passing the site. The CPA time (about ±10 s) and
f0 (about ±15 %) are the reliable parts of a fit.

## What this suggests

- **Events.** Gate on what fits measure well, CPA time and f0, and take speed
  and distance from the real tracks.
- **Seeds.** Two events are often not enough: at McMaster they put the truth
  2nd, and the extras made it 1st. Searching with a few strong events, and
  treating the rest as optional extras on one clock, is the next step.
- **Silence.** As a hard veto between the two largest events it fails on long
  recordings. A softer form would charge a candidate only for flights that
  should clearly have been heard (close, low, loud type, quiet background).
- **Calibration.** "Located" fired falsely three times on development data
  (the railcam twice, Oyster Bay once in a coarse search). The Hamilton
  archive has about 3,000 timestamped recording days at several sites, and
  ADS-B for them is free, which makes it a ready benchmark for calibrating p
  and the hit rate.
