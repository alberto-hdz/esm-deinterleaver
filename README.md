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
- **PRI deinterleaving**: sequential difference histogram (SDIF, Milojevic & Popovic 1992) to find candidate PRIs, then a sequence search that pulls each pulse train out and tolerates dropped pulses

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
- [ ] Jittered and staggered PRIs (PRI transform, stagger detection)
- [ ] Emitter identification against a library
- [ ] Pulse detection from raw I/Q samples
- [ ] Monte Carlo scoring and results

## Known limitations (for now)

- Jittered PRIs aren't separated yet. The fixed sequence-search tolerance can't follow them, so their pulses are reported as one group with no PRI.
- A staggered emitter comes out as several trains at its frame period instead of one emitter.
