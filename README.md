# esm-deinterleaver

Work in progress. I'm building the signal processing chain an ESM (Electronic Support
Measures) receiver or radar warning receiver uses to make sense of the pulses it
intercepts: find the pulses in raw I/Q samples, group them by their parameters,
separate interleaved radars by their pulse timing, and match each one against an
emitter library.

Everything runs on synthetic data from my own simulator, using made-up emitters and
published, open-literature algorithms.

## Pipeline

~~~
raw I/Q -> pulse detector -> PDWs -> DBSCAN clustering -> SDIF PRI analysis -> library ID -> emitters
                             ^
                     (or PDWs straight from the scenario simulator)
~~~

- **Pulse detection**: noise floor from the median of |x|^2 (so the pulses don't bias it), a 9-sample moving average with on/off hysteresis thresholds to find pulses, TOA and PW from interpolated 50% crossings on a lightly smoothed amplitude, and RF from the angle of the summed lag-1 autocorrelation (sum of x[k] times conj(x[k-1]))

- **PDW** (Pulse Descriptor Word): time of arrival, RF, pulse width, amplitude, and angle of arrival for one pulse
- **Clustering**: DBSCAN with an RF-sorted neighbor search, so each lookup only scans pulses within the RF tolerance
- **PRI deinterleaving, pass 1**: sequential difference histogram (SDIF, Milojevic & Popovic 1992) to find candidate PRIs, then a sequence search that pulls each pulse train out and tolerates dropped pulses. Trains that share a PRI inside one cluster are merged into a staggered emitter, and the stagger intervals are recovered from their phases.
- **PRI deinterleaving, pass 2**: for jittered PRIs, a tau-proportional window finds the mean PRI and jitter of the leftover pulses, then the sequence search runs again with a tolerance sized to that jitter
- **Emitter identification**: each emitter is scored against every entry of an emitter library (`config/emitter_library.csv`, all made up) on RF, PW, PRI and PRI type. Values just outside a window still get partial credit, a weak best match is reported as UNKNOWN, and a near-tie with the runner-up is flagged as ambiguous

## Build and run

Needs CMake, a C++17 compiler, GoogleTest, and Python 3 with NumPy.

~~~bash
cmake -S . -B build -G Ninja
cmake --build build
ctest --test-dir build --output-on-failure

python sim/scenario.py --emitters 6 --seed 1 --out data/scn01
./build/esm_cli --pdw data/scn01/pdw.csv --show-truth

# emitters drawn from the library (plus some unknown ones), identified by name
python sim/scenario.py --emitters 8 --seed 3 --library config/emitter_library.csv --out data/lib03
./build/esm_cli --pdw data/lib03/pdw.csv --library config/emitter_library.csv --show-truth

# raw I/Q capture (50 ms at 40 MSPS, 16 MB): detect pulses first, then the same chain
python sim/iq_gen.py --seed 1 --out data/iq01
./build/esm_cli --iq data/iq01/iq.c64 $(cat data/iq01/capture.txt)
~~~

## Status

- [x] Scenario simulator (constant, jittered, and staggered PRIs; pulse loss; spurious pulses)
- [x] Parameter clustering
- [x] SDIF + sequence search for constant-PRI emitters
- [x] Jittered and staggered PRIs
- [x] Emitter identification against a library
- [x] Pulse detection from raw I/Q samples
- [x] Monte Carlo scoring (`sim/score.py`, `sim/score_iq.py`)
- [ ] Result plots and a speed benchmark

## Results so far

Each run is 30 random scenarios, 8 emitters each (2 pairs forced to share RF, AOA and PW), 2% spurious pulses. An emitter counts as found if one output emitter holds at least 80% of its pulses at 90% purity. Reproduce with `python sim/score.py --cli build/esm_cli --runs 30` (add `--drop 0.3` or `--library config/emitter_library.csv`).

| PRI type | found, 10% pulse loss | found, 30% pulse loss |
|---|---|---|
| constant | 120 / 120 | 119 / 120 |
| jitter   | 68 / 76   | 67 / 76   |
| stagger  | 42 / 44   | 43 / 44   |

Identification, with emitters drawn from the library and 25% of them unknown (not in the library), 10% pulse loss:

| | count |
|---|---|
| library emitters given the right name | 206 / 206 |
| unknown emitters reported as UNKNOWN | 31 / 31 |
| wrong names | 0 |

Pulse detection from raw I/Q, 20 captures of 50 ms at 40 MSPS, 4 emitters each at 8 to 25 dB SNR. Reproduce with `python sim/score_iq.py --cli build/esm_cli --runs 20`.

| | result |
|---|---|
| clean pulses detected | 12990 / 12993 (99.98%) |
| false alarms in noise-only time | 0 |
| TOA error, RMS | 7.1 ns (sample period is 25 ns) |
| PW error, RMS | 10.3 ns |
| RF error, RMS | 23.6 kHz |
| pulses overlapping another pulse (merged, not scored above) | 1454 / 14447 (10.1%) |
| emitters found end to end (constant / jitter / stagger) | 37/39, 18/22, 18/19 |

## Known limitations

- Two jittered emitters in the same cluster usually aren't separated; the jitter pass needs most of the leftover gaps to agree on one PRI.
- Two constant emitters with exactly the same PRI in one cluster get merged into a 2-level stagger. From timing alone those look identical.
- A 2-level stagger whose intervals are within a few percent of each other can be reported as jitter.
- The I/Q path has one antenna and no channelizer, so two pulses that overlap in time come out as one merged pulse with the wrong width and frequency. About 10% of pulses overlap in the 4-emitter captures.
