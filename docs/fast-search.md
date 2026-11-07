# Fast search

How the search of this fork differs from dopplertrack 0.3, why, and what it
was measured to cost and gain. The README has the short version.

The work was done in two steps. A research prototype was written first, on
top of an unmodified 0.3 checkout, and measured on the three development
recordings (Oyster Bay, the eagle cam and the Atlanta railcam) on
2026-11-05. It was then ported into the package as the only code path, with
some changes found during the port (an exact closest approach on curved
tracks, a lossless pre-gate, a faster silence filter). Both sets of numbers
are below, labelled. Where they disagree, the package numbers describe the
code in this repository.

Machine for everything on this page: Apple M3 Pro, 12 cores, 18 GB RAM,
11 worker processes.

## What changed

1. **Track store instead of a 1 Hz day index.** Every airborne stretch of
   every flight is kept as knots every 15 s. Each knot has time, position and
   velocity in one Earth-centred frame, as float32. Between two knots the
   track is the cubic Hermite curve through them. Positions and velocities
   come from a 31 s quadratic least-squares fit (Savitzky-Golay) of the 1 Hz
   track, so range rates are not differenced from noisy positions. A 10 s
   light index (one sample per flight per 10 s bucket) sits beside the knots
   for the join and the silence filter. All arrays are `.npy` files opened as
   memory maps.
2. **Exact closest approach on the curve.** `cpa.closest_approach` finds the
   nearest chord among all chords in the window, then runs Newton steps on
   the Hermite arc of every chord that could still be closer (a bound per
   chord rules out the rest). For long windows it starts from the 8 nearest
   chord midpoints and proves that no other chord can win, and falls back to
   all chords when it cannot. On straight, turning and circling test tracks
   it is within 0.3 m of brute force (`closest_approach_dense`).
3. **Log-frequency scorer.** Each event window is resampled once onto a grid
   of query time by log frequency (0.1 % steps), with all harmonics summed.
   An f0 scale is then an index offset, and each placement of a modelled
   curve is one integer gather. 0.3 sampled the linear axis bilinearly for
   every harmonic and every scale.
4. **Two-stage scoring.** Over the seconds that shape a Doppler curve, most
   flights fly nearly straight. So each event gets a table of straight-pass
   scores over CPA distance × speed × time shift, and every candidate gets a
   table lookup. Only the best 20 candidates of the best 500 pairs (by table
   score) are scored exactly on their real tracks.
5. **50 m refinement.** Around the best point of each of the best 50 pairs, a
   50 m lattice (±300 m) is gated the same way and scored exactly.

Removed: the 1 Hz day index and its loader (`DayIndex`, `Segment`,
`build_segments`), scoring every candidate exactly, and the linear-axis
scorer. `dopplertrack store DAY.npz DIR` converts a 0.3 index once.

## Why each step

**Storage.** A 0.3 day index for a busy CONUS day is 6.2 GB and is read into
RAM. On an 18 GB laptop that leaves little room for the worker processes,
and the prototype benchmark ran with heavy swapping (the Oyster Bay lens
stage took 540 s under swap against 110 s without). A store opens in under a
second and only the pages a search touches are read.

**Scoring.** In 0.3, exact scoring of all candidates was most of the time
(16 min of 18 on the upstream author's 16-thread desktop for Oyster Bay).
The log-axis scorer alone cut it 2-3.7× on this Mac. The tables cut it to a
few seconds, because about 1 % of candidates are scored exactly.

**Refinement.** The Oyster Bay peak is sharp: a lattice point 190 m from the
peak scores 24.2 against 29.2 at the peak. With a 200 m lattice, whether the
true pair wins depends on where the lattice happens to fall. The 50 m
refinement removes that luck for about 1-3 s.

## Storage

Store sizes are what `dopplertrack store` writes, light index included,
because a search needs both.

| day | 1 Hz index (0.3) | store: knots | store: light index | store total | ratio |
|---|---|---|---|---|---|
| Oyster Bay, 2026-09-08 | 6.21 GB | 288 MB | 358 MB | 646 MB | 9.6× |
| eagle cam, 2024-12-28 | 1.93 GB | 90 MB | 111 MB | 201 MB | 9.6× |
| railcam, 2026-05-28 | 5.54 GB | 257 MB | 319 MB | 577 MB | 9.6× |
| McMaster Forest, 2026-05-27 | 4.69 GB | 218 MB | 270 MB | 488 MB | 9.6× |
| Safari Road, 2023-07-09 | 1.93 GB | 91 MB | 111 MB | 202 MB | 9.6× |

The knots alone are 21-22× smaller than the index. The light index is
larger than the knots: 14.9 million 10 s samples on the Oyster Bay day, 24
bytes each (position, speed, segment and bucket). The prototype counted
only the knots (328 MB with its own extra arrays) and kept the light index
as a separate 0.84 GB cache; the "19× smaller" figure from that prototype
leaves the cache out.

Conversion from a 0.3 index reads the `.npz` members through memory maps.
It took 17 s (eagle cam), 23 s (Safari Road), 50 s (McMaster), 129 s
(railcam) and 146 s (Oyster Bay), with a peak memory footprint of 1.3-2.8 GB,
on a machine at load averages of 8-18 from other jobs.

Accuracy of the stored curve, measured on the prototype store (same 15 s
knots, same 31 s fit, same knot count): position against the 1 Hz track,
median 2.1-3.4 m and 99th percentile 48-84 m over the three days. The large
errors are in turns. In a standard-rate turn (3°/s at 60 m/s) the 31 s fit
cuts the corner by about 8 m and the speed by about 4 m/s.

## Package measurements

Measured on 2026-11-06 with `dopplertrack search ... --workers 11` (the
commands are in [benchmarks/](../benchmarks/)), one search at a time, but
**not on an idle machine**: other jobs held the load average at 19-60 when
each search started (two other CPU-bound processes ran throughout). Wall
times are therefore upper bounds, and CPU time (user, all processes) is the
steadier figure. Peak memory is the main process's peak footprint. These
runs predate the 0.25 s shift grid and the lens prefilter described under
[Validation against the prototype](#validation-against-the-prototype-2026-11-06),
which has runs of the current code; the held-out rows there were rerun too
(McMaster 89 s and Safari Road 152 s with `--no-silence`, at load averages
of 63-68, same winners).

| recording | joined pairs | pairs with candidates | candidates | join | lens | tables + lookup | exact | refine | total (wall) | CPU | peak memory |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Oyster Bay | 60,802 | 9,448 | 2.44 M | 8.0 s | 72.1 s | 2.4 s | 1.0 s | 2.1 s | 86 s | 664 s | 1.7 GB |
| eagle cam | 96,814 | 12,544 | 1.04 M | 4.2 s | 45.7 s | 2.4 s | 1.1 s | 1.0 s | 54 s | 458 s | 1.4 GB |
| railcam (`--no-silence`) | 179,222 | 32,738 | 1.49 M | 9.0 s | 61.8 s | 2.6 s | 1.1 s | 0.9 s | 75 s | 577 s | 1.6 GB |
| McMaster (`--no-silence`) | 92,922 | 19,190 | 1.84 M | 6.9 s | 46.9 s | 2.8 s | 1.2 s | 1.0 s | 59 s | 456 s | 1.7 GB |
| McMaster, silence on | 92,922 | 481 | 31,631 | 9.3 s | 142.3 s | 2.2 s | 1.0 s | 3.7 s | 159 s | 1,181 s | 1.5 GB |
| Safari Road (`--no-silence`) | 174,948 | 33,055 | 0.99 M | 5.1 s | 104.1 s | 2.4 s | 1.1 s | 0.7 s | 113 s | 1,001 s | 2.4 GB |
| Safari Road, silence on | 174,948 | 985 | 21,309 | 5.1 s | 257.5 s | 2.3 s | 1.0 s | 1.3 s | 267 s | 1,292 s | 2.4 GB |

Totals exclude the extras (under 1 s) and loading the store and spectrogram
(0.2-2.4 s). The lens stage is 80-96 % of every search.

The silence filter is the slow part of the lens when there is a long
stretch between the two events (5.3 h at McMaster, 47 min at Safari Road):
almost every point has other flights to check. Skipping points already silenced cut the McMaster
search from 2,046 s to 1,181 s of CPU with identical results; Safari Road's
wall time rose in the rerun (159 s to 267 s) while its CPU time fell (1,422 s
to 1,292 s), a sign of how much the machine's load moved the wall times.

Two of the same searches with upstream 0.3 (its 1 Hz day index, every
candidate scored exactly with its linear-axis scorer), same machine, same
events, 11 workers, run right after the fork's:

| | upstream 0.3 | this fork | ratio |
|---|---|---|---|
| **eagle cam**, wall time | 309 s (index load 2 s, join 4 s, lens 35 s, scoring 258 s) | 58 s | 5.4× |
| eagle cam, CPU time | 2,654 s | 458 s | 5.8× |
| eagle cam, peak memory | 3.6 GB | 1.4 GB | |
| **Oyster Bay**, wall time | 2,043 s (index load 5 s, join 8 s, lens 94 s, scoring 1,865 s) | 88 s | 23× |
| Oyster Bay, CPU time | 12,086 s | 664 s | 18× |
| Oyster Bay, peak memory | 8.5 GB | 1.7 GB | |

The fork's times include the extras; upstream's include the extras and its
index load. Two things make the Oyster Bay ratio larger than the eagle
cam's. The 6.2 GB index pushed the 18 GB machine into swap (10-12 GB of swap
in use), and upstream's scorer slowed from 0.2 ms to 0.7 ms per candidate.
And the test suite ran during 7 of upstream's 34 minutes. On a machine with RAM to
spare, expect a ratio nearer the eagle cam's. The prototype had measured
the full log-axis search at 418 s for Oyster Bay on a quiet machine, 7× its
own 60 s.

Accuracy of the same two runs (the fork with 0.5 s shift steps; the current
numbers are in the next table; the fork's p-values moved by up to 0.32):

| | upstream 0.3 | this fork |
|---|---|---|
| Oyster Bay, two events | true pair 1 of 9,669; 32 m; 6.2 % ahead; ambiguous, p = 0.23 | 1 of 498; 26 m; 8.1 % ahead; located, p = 0.09 |
| Oyster Bay, with the extra | 1; 32 m; 24.3 % ahead; located, p = 0.0001 | 1; 96 m; 11.3 % ahead; located, p = 0.05 |
| eagle cam, two events | 13 of 13,302; 972 m; winner 203 km off; no match, p = 0.52 | 13 of 500; 800 m; winner 203 km off; no match, p = 0.94 |
| eagle cam, with extras | 1; 627 m; 13.5 % ahead; located, p = 0.05 | 1; 741 m; 4.9 % ahead; ambiguous, p = 0.42 |

Same winners and the same ranks of the true pair. The fork's margins with
extras are smaller and its verdicts weaker: its null is the 500
shortlisted pairs, all of them good fits, where upstream's is every pair.

Accuracy of the package on the development recordings, from the
back-to-back runs below (rank of the true pair among the best pair per
flight pair; distance of its best point from the answer):

| | two events | with extra flyovers |
|---|---|---|
| Oyster Bay | 1 of 499; 26 m; 8.1 % ahead; located, p = 0.06 | 1 of 199; 96 m; 12.4 % ahead; located, p = 0.04 |
| eagle cam | 15 of 500 (797 m); the winner is 203 km off; no match, p = 0.62 | 1 of 200; 741 m; 4.9 % ahead; ambiguous, p = 0.43 |
| railcam | 4 of 500 (95 m); the winner is 692 km off; ambiguous, p = 0.31 | — |

## Validation against the prototype (2026-11-06)

Each recording was searched with this package and, right after, with the
research prototype, on the same machine with the same events and 11
workers. Stores were converted from the 0.3 indexes with `dopplertrack
store`; the prototype used its own stores. Load averages at the start of
each run were 7-22 (another agent's CPU-bound jobs ran throughout). The
prototype's times leave out the second exact stage it ran and discarded
(K = 1000, N = 50; 6.2-6.9 s).

| | Oyster Bay | eagle cam | railcam (`--no-silence`) |
|---|---|---|---|
| package: join / lens / tables + lookup / exact / refine | 8.1 / 64.0 / 4.1 / 1.6 / 2.5 s | 4.0 / 41.9 / 4.0 / 1.8 / 1.5 s | 9.2 / 53.9 / 4.3 / 1.8 / 1.3 s |
| package: stages / wall / CPU | 80 s / 82 s / 637 s | 53 s / 56 s / 462 s | 70 s / 73 s / 571 s |
| prototype: join / lens / tables + lookup / exact / refine | 7.8 / 44.4 / 2.5 / 1.7 / 1.1 s | 4.2 / 24.1 / 2.7 / 1.7 / 1.1 s | 8.6 / 28.8 / 2.6 / 1.6 / 1.0 s |
| prototype: stages (incl. load) / CPU | 60 s / 496 s | 35 s / 335 s | 45 s / 359 s |
| peak memory footprint, package / prototype | 2.08 / 2.79 GB | 1.44 / 1.93 GB | 1.73 / 2.73 GB |
| store on disk, package / prototype (with its light cache) | 646 MB / 1,463 MB | 201 MB / 530 MB | 577 MB / 1,252 MB |
| true pair, package / prototype | 1st / 1st | 15th / 13th; 1st / 1st with extras | 4th / 3rd |

The package's lens stage is 1.4-1.9× the prototype's. Its closest
approaches are exact on the Hermite curve (chords, then Newton steps on
every arc that could be nearer), where the prototype took the nearest knot
and the two chords beside it. On the gate calls of 1,500 Oyster Bay pairs
the exact closest approach took 3.0 s against 0.95 s for the chord step
alone (flight 1, every exact gate call). On 3,000-pair samples run in one
process, the package's lens used 1.4× (Oyster Bay), 1.7× (eagle cam) and 1.9× (railcam) the prototype's CPU
time; the prototype found 0.7 % and 1.9 % more candidates on the first two
(its silence filter, on chords, silences fewer points) and 2.3 % fewer on
the railcam (its pre-gate). Three things were tried to close the gap:

- a chord-distance prefilter before the exact coarse gate (kept: identical
  candidates on 4,000-pair samples of each day, 7-17 % less lens CPU);
- computing each chain's chords once (kept: 13 % on a railcam sample);
- judging the coarse gate on the chords, or on bounds from the chords, with
  the curve only where two chords contend (dropped: the first lost 0.02 %
  of the candidates, the second let 60 % more fine points through and was
  no faster).

The runs in the table are from the commit before the first two of these.
The same six runs repeated with the final code, at load averages of 42-70,
gave identical ranks, scores and verdicts, with lens ratios of 1.5-2.2 and
CPU ratios of 1.2-1.6; their wall times were 1.5-1.6 times those above.

**The railcam rank.** With 0.5 s time-shift steps the package ranked the
railcam's true pair 53rd, the prototype 3rd. At the same lattice point the
package scored it 11.17 and the prototype 11.47. The difference was the
heard CPA time: the package's exact closest approach put it 0.10 s later
than the prototype's chords. Event 1's best shift then moved from 0.5 s to
1.0 s, and event 2, which takes event 1's shift, lost 0.24 of its score.
That left the pair 53rd after the exact stage, below the 50 that get
refined, where the prototype had it 26th. With 0.25 s steps (now the
default) it is 33rd after the exact stage and 4th after refinement. Rescoring
the package's shortlist offline, choosing the shift by the sum of both
events instead also helped (37th), but changes the method; it was not
adopted. The finer grid costs about 2 s per search (tables and exact
scores) and moved the eagle cam's two-event rank from 13th to 15th.

**Store size including caches.** A store holds everything a search reads
from the ADS-B day; there is no other cache. Terrain is separate: decoded
zoom-9 tiles are kept under `~/.cache/dopplertrack/terrarium/9` (PNG plus a
float32 `.npy` per tile) and shared by every day; after all the searches on
this page it held 6,801 tiles, 2.35 GB (1.81 GB of it `.npy`). `refine`
uses zoom-12 tiles in their own cache.

**Store conversion** (`dopplertrack store`, 2026-11-06, load 9-22): Oyster Bay
62 s, eagle cam 20 s, railcam 54 s, with peak memory footprints of 2.77, 1.25
and 2.57 GB.

## Prototype measurements (2026-11-05)

These came from the research prototype on a quiet machine, before the
port. "Original" is 0.3's linear-axis scorer; its times on all candidates
were **estimated** from a 4,000-row sample, not measured, and were taken
while the machine was swapping. "Log-axis full" scores every candidate with
the log-axis scorer, and was measured.

Time per search, s:

| stage | eagle cam: full / fast | railcam: full / fast | Oyster Bay: full / fast |
|---|---|---|---|
| index or store load | 0.6 / 0.7 | 3.2 / 2.1 | 3.9 / 0.6 |
| join and lens | 32 / 26 | 44 / 35 | 91 / 54 |
| score tables | — / 2.4 | — / 3.1 | — / 1.9 |
| scoring: original (est.) / log-axis full / tables, exact and refine | ~300 / 152 / 3 | ~600 / 210 / 3 | ~1440 / 323 / 3.5 |
| total | ~333 (log-axis 185) → 32 | ~650 (log-axis 257) → 43 | ~1535 (log-axis 418) → 60 |
| speed-up vs original (est.) / vs log-axis full | 10× / 5.8× | 15× / 6× | 26× / 7× |

The 10-26× figures are against the estimated original scorer. The measured
speed-up is the second one, 5.8-7× against the full log-axis search.

Rank of the true pair:

| | eagle cam | railcam | Oyster Bay |
|---|---|---|---|
| log-axis, all candidates scored exactly | 13 | 53; a wrong pair "located", p = 0.04 | 1; 7.0 % ahead; p = 0.20; 32 m |
| tables, exact best 20 × 500, no refinement | 7 | 26 | 22 |
| fast, two events (with refinement) | 13 | 3; a wrong pair wins, "ambiguous", p = 0.25 | 1; 8.0 % ahead; p = 0.08; about 26 m |
| fast with extra flyovers | 1; 5.0 % ahead; p = 0.32 | — | 1; 11.2 % ahead; p = 0.04 |

The refinement moved the Oyster Bay truth from 22nd to 1st and the railcam
truth from 26th to 3rd, but the eagle-cam truth from 7th to 13th. On the
railcam the winner is still a wrong pair.

Candidates scored exactly: about 10,000 of 1.0-2.5 million (0.4-1.0 %), plus
6,300-7,700 refinement points.

## What did not work

These were tried in the prototype and dropped. Numbers from 2026-11-05.

- **Table gates in the lens stage.** Dropping candidates whose table score was
  below half the table maximum cut candidates 5-12× on Oyster Bay and the
  railcam, saved only 13 s, and lost the Oyster Bay truth. The high-scoring
  region of a table is a diagonal ridge (distance roughly proportional to
  speed); without refinement, the best points of a sharp peak fall outside
  it.
- **Coarse-then-fine by random subsampling.** Scoring a random 1/8 of
  candidates and then the top 300 pairs in full lost the truth on both
  sharp-peaked recordings (Oyster Bay rank 479, railcam 663).
- **A gate on aircraft type from the rest frequency.** A third of flights have
  no type, so it removes at most about 2× in busy airspace. It dropped the
  railcam truth, because `detect` had fitted 1.5× the blade tone.
- **Jet loudness as a location score.** The 400-1500 Hz level over the hour
  correlates 0.36 with the ADS-B prediction at the Oyster Bay answer, against
  0.14 for the runner-up. It was noise at the eagle cam (0.03), and a wrong
  pair in Nebraska also reached 0.32. Usable as a check on a shortlist, not
  alone. See [other-signals.md](other-signals.md).

## Caveats

- The settings (500 pairs, 20 points, 50 refined pairs) were chosen after
  seeing Oyster Bay fail without refinement, on the same three recordings.
  Only the railcam result is partly independent of that choice. The
  held-out recordings ([heldout.md](heldout.md)) are the independent test.
- The verdict's wrong pairs are the shortlisted pairs only, some refined and
  some not. Its p-values are not calibrated for the two-stage search; the CLI
  prints a note saying so.
- The lens stage (join, lattices, closest approaches, gates, silence filter)
  is now 77-92 % of the time. Most joined pairs produce no candidates (9,448
  of 60,802 on the Oyster Bay day do). Rejecting those earlier is the next
  speed step.
