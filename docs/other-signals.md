# Other signals in the recordings

What else the Oyster Bay and eagle-cam recordings carry besides propeller
Doppler curves, tested on 2026-11-05 with throwaway research scripts that
are not part of this package. Ground truth is the location in
[experiments.md](experiments.md). Oyster Bay starts at 20:17:29.8 UTC on
2026-09-08. The eagle cam's start, 14:21:09 UTC on 2024-12-28, was recovered
by matching ADS-B passes to `detect` events. ADS-B came from adsb.lol for
both days, within about 50 km of each site, jets included.

One or two recordings per idea, ground truth used throughout, no blind
runs. Treat these as leads, not results.

| idea | result |
|---|---|
| Jets, broadband loudness | **Works at Oyster Bay.** The 400-1500 Hz level over the hour correlates 0.40 with the ADS-B prediction (sum of 1/r² with air absorption over every aircraft) at the true site. At sites moved away: 0.36 at 1 km, 0.25 at 2 km, 0.13 at 3 km, about 0 beyond 5 km. With the place known, sliding over the whole day puts the start at the right time (10 s grid), but the runner-up is close (0.41 against 0.36, z = 4.5). **Fails on the eagle cam** (0.03): wind and birds, few jets. |
| Jets, Doppler tones above 450 Hz | Tonal curves at 0.4-2.5 kHz are visible on some low jets (a Boeing 767, 2000 to 1400 Hz). Locking on them is unreliable: 2 of 16 passes clearly beat the null. |
| Aircraft type from the rest frequency | Consistent with the type data. Three PC-12s sit at 111.0, 111.7 and 112.4 Hz (4 blades at about 1665-1686 rpm, within 1.3 %). A King Air 200 is at 112.2 Hz and a Twin Comanche at 83.9 Hz (2 blades, about 2517 rpm). A King Air 350 at 125.7 Hz implies 4 blades at 1886 rpm (above the usual 1700) or 5 blades at 1508 rpm; the blade count is unconfirmed. Turboprops are tight; pistons vary with power setting. No shaft-rate lines were visible. |
| Loudness and harmonic tilt as range | **Not usable as built.** The envelopes are limited by SNR (only about ±50 s above the noise) and fall far faster than 1/r²: excess attenuation at low elevation dominates. Range fits run to the edge of the allowed range. |
| Ground-reflection fringes | **Inconclusive.** Fitted microphone heights disagree across passes (0.8-2.0 m), and contrast is only about 2× a shuffled-geometry null. Needs a soft-ground model. |
| Fore/aft loudness asymmetry | Approach is louder on the PC-12 (+4 to +5 dB on harmonics 1-3) and the King Air 350 (+2.6 dB on harmonic 1, peak about 25 s before CPA); the Twin Comanche shows about 0 dB. Too few events to call it a feature. |
| Stereo bearing | **Dead on both.** Oyster Bay has identical channels. The eagle cam's channels differ by 10 dB, but their coherent part sits at zero lag through every flyover: one source mixed differently, not two microphones. |
| Power changes against climb or descent | Not tested; the passes found were in level flight. |
| Atmosphere (speed of sound, wind) | Not tested. A 2 °C error changes c by 0.35 %, which the f0 fit absorbs almost entirely. Too weak per pass. |
| Detection range | Propeller tones locked at 3.6, 6.1 and 8.2 km at Oyster Bay. On the eagle cam only 13 of 101 passes within 6 km beat the null. The lock statistic is crude, so this is a lower bound. |
| One clock for all events | Oyster Bay's propeller passes agree within about 1 s (+0.75, -0.5 to +0.25, +0.75 s); jet loudness says -2 s. Eagle-cam per-pass shifts scatter ±4 s. |
| Mains hum (ENF) | **Oyster Bay has it**, though recorded outdoors: the 120, 240 and 360 Hz lines wander together (360 Hz at +17 dB in 93 % of 16 s frames, ±0.03 Hz; the 240 and 360 Hz tracks correlate 0.54), so the source is locked to the grid. Not yet matched to a grid-frequency reference. **The eagle cam has none** on either channel. |

The jet-loudness check is the only one tried against a search result: on
the 200 best Oyster Bay pairs it ranks the true pair 1st (0.36, against 0.14
for a runner-up pair over Texas). It was noise at the eagle cam and weak at
the railcam, and a wrong pair in Nebraska also reached 0.32. It is not in
the package.
