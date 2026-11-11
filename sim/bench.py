"""Time the C++ pipeline.

1. Real-time check: a run of 200 ms windows with 16 emitters each, the way a
   receiver would hand the processor one buffer at a time. Reports pulses per
   second and how much faster than real time the processing runs.
2. Scaling: one long capture at a time (up to 10 s / ~340k pulses) to see how
   clustering and PRI analysis grow with the number of pulses.

Timings come from esm_cli's own stage timers, so file loading isn't counted.

usage: python sim/bench.py --cli build/esm_cli
"""
import argparse
import re
import subprocess
import sys
import tempfile
from pathlib import Path

SIM = Path(__file__).resolve().parent / "scenario.py"
TIMING = re.compile(r"pipeline .* ([\d.]+) ms \(clustering ([\d.]+), PRI ([\d.]+)\)")


def run(cli, out, duration_ms, seed, emitters=16):
    subprocess.run([sys.executable, str(SIM), "--emitters", str(emitters), "--same-rf-pairs", "4",
                    "--duration-ms", str(duration_ms), "--seed", str(seed), "--out", str(out)],
                   check=True, capture_output=True)
    pulses = sum(1 for _ in open(out / "pdw.csv")) - 1
    text = subprocess.run([cli, "--pdw", str(out / "pdw.csv")], check=True, capture_output=True,
                          text=True).stdout
    total, clus, pri = (float(v) for v in TIMING.search(text).groups())
    return pulses, total, clus, pri


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--cli", default="build/esm_cli")
    p.add_argument("--windows", type=int, default=20)
    a = p.parse_args()

    with tempfile.TemporaryDirectory() as work:
        print(f"real-time check: {a.windows} windows of 200 ms, 16 emitters each")
        pulses = total = 0.0
        for seed in range(a.windows):
            n, t, _, _ = run(a.cli, Path(work) / f"w{seed}", 200, seed)
            pulses += n
            total += t
        data_ms = 200.0 * a.windows
        print(f"  {int(pulses)} pulses processed in {total:.0f} ms "
              f"-> {pulses / total * 1000:,.0f} pulses/s, {data_ms / total:.1f}x faster than real time")

        print("scaling (one capture each, 16 emitters):")
        print(f"  {'duration':>9} {'pulses':>8} {'clustering ms':>14} {'PRI ms':>9}")
        for dur in (200, 1000, 2000, 10000):
            n, _, clus, pri = run(a.cli, Path(work) / f"s{dur}", dur, 7)
            print(f"  {dur:>7} ms {n:>8} {clus:>14.1f} {pri:>9.1f}")


if __name__ == "__main__":
    main()
