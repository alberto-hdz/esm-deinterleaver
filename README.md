# esm-deinterleaver

[![ci](https://github.com/alberto-hdz/esm-deinterleaver/actions/workflows/ci.yml/badge.svg)](https://github.com/alberto-hdz/esm-deinterleaver/actions/workflows/ci.yml)

The signal processing chain an ESM (Electronic Support Measures) receiver or radar
warning receiver uses to make sense of the pulses it intercepts, written in C++17:
find the pulses in raw I/Q samples, group them by their parameters, separate
interleaved radars by their pulse timing, and match each one against an emitter
library.

Everything runs on synthetic data from my own Python simulators, using made-up
emitters and published, open-literature algorithms.

## Pipeline

~~~
raw I/Q -> pulse detector -> PDWs -> DBSCAN clustering -> SDIF PRI analysis -> library ID -> emitters
                             ^
                     (or PDWs straight from the scenario simulator)
~~~

- **Pulse detection**: noise floor from the median of |x|^2 (so the pulses don't bias it), a 9-sample moving average with on/off hysteresis thresholds to find pulses, TOA and PW from interpolated 50% crossings on a lightly smoothed amplitude, and RF from the angle of the summed lag-1 autocorrelation (sum of x[k] times conj(x[k-1]))
- **PDW** (Pulse Descriptor Word): time of arrival, RF, pulse width, amplitude, and angle of arrival for one pulse
- **Clustering**: DBSCAN on RF, AOA and PW. The pulses are sorted by RF, the core test stops as soon as it has enough neighbors, and a skip list over the RF order means growing a cluster only visits pulses that can still join it
- **PRI deinterleaving, pass 1**: sequential difference histogram (SDIF, Milojevic & Popovic 1992) to find candidate PRIs, then a sequence search that pulls each pulse train out and tolerates dropped pulses. A least-squares fit through each train gives the exact PRI. Trains that share a PRI inside one cluster are merged into a staggered emitter, and the stagger intervals are recovered from their phases. Trains whose missing pulses repeat on a schedule are rejected, since real pulse loss is random
- **PRI deinterleaving, pass 2**: for jittered PRIs, a tau-proportional window finds the mean PRI and jitter of the leftover pulses, then the sequence search runs again with a tolerance sized to that jitter
- **Emitter identification**: each emitter is scored against every entry of an emitter library (`config/emitter_library.csv`, all made up) on RF, PW, PRI and PRI type. Values just outside a window still get partial credit, a weak best match is reported as UNKNOWN, and a near-tie with the runner-up is flagged as ambiguous

## Build and run

Needs CMake, a C++17 compiler, and Python 3 with NumPy (and matplotlib for the plots).
GoogleTest is used if it's installed, otherwise CMake downloads it.

~~~bash
cmake -S . -B build -G Ninja
cmake --build build
ctest --test-dir build --output-on-failure
python -m pytest -q sim/tests

python sim/scenario.py --emitters 6 --seed 1 --out data/scn01
./build/esm_cli --pdw data/scn01/pdw.csv --show-truth

# emitters drawn from the library (plus some unknown ones), identified by name
python sim/scenario.py --emitters 8 --seed 3 --library config/emitter_library.csv --out data/lib03
./build/esm_cli --pdw data/lib03/pdw.csv --library config/emitter_library.csv --show-truth

# raw I/Q capture (50 ms at 40 MSPS, 16 MB): detect pulses first, then the same chain
python sim/iq_gen.py --seed 1 --out data/iq01
./build/esm_cli --iq data/iq01/iq.c64 $(cat data/iq01/capture.txt)
~~~

## Results

An emitter counts as found if one output emitter holds at least 80% of its pulses at
90% purity. Each scenario has 2 pairs of emitters forced to share RF, AOA and PW, so
only the PRI analysis can separate them, plus 2% spurious pulses.

![Emitters found vs pulse loss](docs/found_vs_loss.png)

![Emitters found vs number of emitters](docs/found_vs_emitters.png)

30 scenarios of 8 emitters (`python sim/score.py --cli build/esm_cli --runs 30`, add `--drop 0.3`):

| PRI type | found, 10% pulse loss | found, 30% pulse loss |
|---|---|---|
| constant | 120 / 120 | 119 / 120 |
| jitter   | 68 / 76   | 67 / 76   |
| stagger  | 42 / 44   | 43 / 44   |

Identification, with emitters drawn from the library and 25% of them unknown, 10% pulse loss (add `--library config/emitter_library.csv`):

| | count |
|---|---|
| library emitters given the right name | 206 / 206 |
| unknown emitters reported as UNKNOWN | 31 / 31 |
| wrong names | 0 |

### Raw I/Q path

![Pulse detection and TOA error vs SNR](docs/iq_vs_snr.png)

20 captures of 50 ms at 40 MSPS, 4 emitters each at 8 to 25 dB SNR (`python sim/score_iq.py --cli build/esm_cli --runs 20`):

| | result |
|---|---|
| clean pulses detected | 12990 / 12993 (99.98%) |
| false alarms in noise-only time | 0 |
| TOA error, RMS | 7.1 ns (sample period is 25 ns) |
| PW error, RMS | 10.3 ns |
| RF error, RMS | 23.6 kHz |
| pulses overlapping another pulse (merged, not scored above) | 1454 / 14447 (10.1%) |
| emitters found end to end (constant / jitter / stagger) | 37/39, 18/22, 18/19 |

The plots come from `python sim/plots.py --cli build/esm_cli --out docs`.

## Performance

`python sim/bench.py --cli build/esm_cli`, single thread, Release build. Processing
200 ms windows with 16 emitters each, the pipeline handles about 380,000 pulses per
second, roughly 9x faster than real time.

Clustering before and after the rewrite (same output on 80 test scenarios):

| capture | pulses | clustering before | clustering after |
|---|---|---|---|
| 0.2 s  | 6.8k  | 93 ms     | 1.2 ms |
| 1 s    | 34k   | 2,380 ms  | 7.9 ms |
| 2 s    | 67k   | 13,974 ms | 22 ms  |
| 10 s   | 337k  | out of memory | 254 ms |

The first version was O(n^2): every pulse in a big cluster re-scanned its whole RF
window, and a pulse could be queued once per neighbor, which ran out of memory at
337k pulses. These numbers are from a Linux VM; `bench.py` prints the ones for your
own machine.

## Testing

- 56 GoogleTest unit tests (PDW I/O, clustering, PRI analysis, identification, pulse detection) and 11 pytest tests for the simulators
- Several tests come from bugs I hit along the way, like a constant emitter turning into a fake stagger at 30% pulse loss, weak pulses splitting in two, and a PRI shorter than the search tolerance looping until it ran out of memory
- AddressSanitizer / UBSan build:
  `cmake -S . -B build-asan -DCMAKE_BUILD_TYPE=Debug -DCMAKE_CXX_FLAGS="-fsanitize=address,undefined -fno-omit-frame-pointer"`
- clang-tidy clean with the checks in `.clang-tidy`: `clang-tidy -p build src/*.cpp apps/*.cpp`
- GitHub Actions builds and runs everything on every push, including a smoke test of all three input modes

## Known limitations

- Two jittered emitters in the same cluster usually aren't separated; the jitter pass needs most of the leftover gaps to agree on one PRI.
- Two constant emitters with exactly the same PRI in one cluster get merged into a 2-level stagger. From timing alone those look identical.
- A 2-level stagger whose intervals are within a few percent of each other can be reported as jitter.
- The I/Q path has one antenna and no channelizer, so two pulses that overlap in time come out as one merged pulse with the wrong width and frequency. About 10% of pulses overlap in the 4-emitter captures.
- PRI analysis time grows faster than linearly with very long captures (7.7 s for 337k pulses), because more PRI candidates clear the threshold. In 200 ms windows it isn't a problem.
