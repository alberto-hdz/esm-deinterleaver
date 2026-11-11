# esm-deinterleaver

Work in progress. I'm building the signal processing chain an ESM (Electronic Support
Measures) receiver or radar warning receiver uses to make sense of the pulses it
intercepts: group the pulses by their parameters, separate interleaved radars by
their pulse timing, and (next) match each one against an emitter library.

Everything runs on synthetic data from my own simulator, using made-up emitters and
published, open-literature algorithms.

## Pipeline

~~~
pulses (PDWs) -> DBSCAN clustering on RF / AOA / PW -> SDIF PRI analysis + sequence search -> emitters
~~~

- **PDW** (Pulse Descriptor Word): time of arrival, RF, pulse width, amplitude, and angle of arrival for one pulse
- **Clustering**: DBSCAN with an RF-sorted neighbor search, so each lookup only scans pulses within the RF tolerance
- **PRI deinterleaving, pass 1**: sequential difference histogram (SDIF, Milojevic & Popovic 1992) to find candidate PRIs, then a sequence search that pulls each pulse train out and tolerates dropped pulses. Trains that share a PRI inside one cluster are merged into a staggered emitter, and the stagger intervals are recovered from their phases.
- **PRI deinterleaving, pass 2**: for jittered PRIs, a tau-proportional window finds the mean PRI and jitter of the leftover pulses, then the sequence search runs again with a tolerance sized to that jitter

## Build and run

Needs CMake, a C++17 compiler, GoogleTest, and Python 3 with NumPy.

~~~bash
cmake -S . -B build -G Ninja
cmake --build build
ctest --test-dir build --output-on-failure

python sim/scenario.py --emitters 6 --seed 1 --out data/scn01
./build/esm_cli --pdw data/scn01/pdw.csv --show-truth
~~~

## Status

- [x] Scenario simulator (constant, jittered, and staggered PRIs; pulse loss; spurious pulses)
- [x] Parameter clustering
- [x] SDIF + sequence search for constant-PRI emitters
- [x] Jittered and staggered PRIs
- [ ] Emitter identification against a library
- [ ] Pulse detection from raw I/Q samples
- [ ] Monte Carlo scoring and results

## Results so far

30 random scenarios, 8 emitters each (2 pairs forced to share RF, AOA and PW), 10% pulse loss, 2% spurious pulses. An emitter counts as found if one output emitter holds at least 80% of its pulses at 90% purity. Reproduce with `python sim/score.py --cli build/esm_cli --runs 30`.

| PRI type | found |
|---|---|
| constant | 120 / 120 |
| jitter   | 68 / 76 |
| stagger  | 42 / 44 |

## Known limitations

- Two jittered emitters in the same cluster usually aren't separated; the jitter pass needs most of the leftover gaps to agree on one PRI.
- Two constant emitters with exactly the same PRI in one cluster get merged into a 2-level stagger. From timing alone those look identical.
- A 2-level stagger whose intervals are within a few percent of each other can be reported as jitter.
