"""Score esm_cli against the simulator's ground truth over many random scenarios.

For each seed I generate a scenario, run esm_cli on it, and check every true
emitter: it counts as found if one output emitter holds at least 80% of its
pulses (recall) and at least 90% of that output emitter's pulses are its own
(purity). Results are grouped by PRI type.

usage: python sim/score.py --cli build/esm_cli --runs 30
"""
import argparse
import collections
import csv
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

SIM = Path(__file__).resolve().parent / "scenario.py"


def score_one(cli, seed, emitters, pairs, drop, spurious, workdir):
    out = Path(workdir) / f"s{seed}"
    subprocess.run([sys.executable, str(SIM), "--emitters", str(emitters), "--same-rf-pairs",
                    str(pairs), "--drop", str(drop), "--spurious", str(spurious),
                    "--seed", str(seed), "--out", str(out)], check=True, capture_output=True)
    subprocess.run([cli, "--pdw", str(out / "pdw.csv"), "--labels", str(out / "labels.csv")],
                   check=True, capture_output=True)

    truth = list(csv.DictReader(open(out / "truth.csv")))
    tid = np.array([int(r["truth_id"]) for r in csv.DictReader(open(out / "pdw.csv"))])
    lab = np.array([int(r["emitter"]) for r in csv.DictReader(open(out / "labels.csv"))])

    found = collections.Counter()
    total = collections.Counter()
    for t in truth:
        mine = tid == int(t["id"])
        total[t["pri_type"]] += 1
        labels, counts = np.unique(lab[mine], return_counts=True)
        for L, cnt in zip(labels, counts):
            if L < 0:
                continue
            recall = cnt / mine.sum()
            purity = cnt / (lab == L).sum()
            if recall >= 0.8 and purity >= 0.9:
                found[t["pri_type"]] += 1
                break
    n_out = len(set(lab.tolist()) - {-1})
    return found, total, n_out, len(truth)


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--cli", default="build/esm_cli")
    p.add_argument("--runs", type=int, default=30)
    p.add_argument("--emitters", type=int, default=8)
    p.add_argument("--same-rf-pairs", type=int, default=2)
    p.add_argument("--drop", type=float, default=0.10)
    p.add_argument("--spurious", type=float, default=0.02)
    a = p.parse_args()

    found = collections.Counter()
    total = collections.Counter()
    n_out = n_true = 0
    with tempfile.TemporaryDirectory() as work:
        for seed in range(a.runs):
            f, t, o, n = score_one(a.cli, seed, a.emitters, a.same_rf_pairs, a.drop, a.spurious, work)
            found.update(f)
            total.update(t)
            n_out += o
            n_true += n

    print(f"{a.runs} scenarios, {a.emitters} emitters each, drop {a.drop:.0%}, spurious {a.spurious:.0%}")
    for kind in ("constant", "jitter", "stagger"):
        print(f"  {kind:<9} {found[kind]:4d} / {total[kind]:<4d} found")
    print(f"  emitters reported: {n_out} (true: {n_true})")


if __name__ == "__main__":
    main()
