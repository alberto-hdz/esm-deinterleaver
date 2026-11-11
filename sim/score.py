"""Score esm_cli against the simulator's ground truth over many random scenarios.

For each seed I generate a scenario, run esm_cli on it, and check every true
emitter: it counts as found if one output emitter holds at least 80% of its
pulses (recall) and at least 90% of that output emitter's pulses are its own
(purity). Results are grouped by PRI type.

With --library, emitters come from the library (plus some unknown ones) and I
also check the name esm_cli gave each found emitter against the truth.

usage: python sim/score.py --cli build/esm_cli --runs 30
       python sim/score.py --cli build/esm_cli --runs 30 --library config/emitter_library.csv
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


def score_one(cli, seed, emitters, pairs, drop, spurious, workdir, library, unknown_frac):
    out = Path(workdir) / f"s{seed}"
    sim_cmd = [sys.executable, str(SIM), "--emitters", str(emitters), "--same-rf-pairs",
               str(pairs), "--drop", str(drop), "--spurious", str(spurious),
               "--seed", str(seed), "--out", str(out)]
    cli_cmd = [cli, "--pdw", str(out / "pdw.csv"), "--labels", str(out / "labels.csv"),
               "--emitters", str(out / "emitters.csv")]
    if library:
        sim_cmd += ["--library", library, "--unknown-frac", str(unknown_frac)]
        cli_cmd += ["--library", library]
    subprocess.run(sim_cmd, check=True, capture_output=True)
    subprocess.run(cli_cmd, check=True, capture_output=True)

    truth = list(csv.DictReader(open(out / "truth.csv")))
    tid = np.array([int(r["truth_id"]) for r in csv.DictReader(open(out / "pdw.csv"))])
    lab = np.array([int(r["emitter"]) for r in csv.DictReader(open(out / "labels.csv"))])
    out_name = {int(r["emitter"]): r["name"] for r in csv.DictReader(open(out / "emitters.csv"))}

    found = collections.Counter()
    total = collections.Counter()
    ident = collections.Counter()  # how the names came out, for emitters that were found
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
                got, want = out_name[int(L)], t["name"]
                if want == "UNKNOWN":
                    ident["unknown, called UNKNOWN" if got == "UNKNOWN" else "unknown, given a name"] += 1
                elif got == want:
                    ident["known, right name"] += 1
                elif got == "UNKNOWN":
                    ident["known, called UNKNOWN"] += 1
                else:
                    ident["known, wrong name"] += 1
                break
    n_out = len(set(lab.tolist()) - {-1})
    return found, total, ident, n_out, len(truth)


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--cli", default="build/esm_cli")
    p.add_argument("--runs", type=int, default=30)
    p.add_argument("--emitters", type=int, default=8)
    p.add_argument("--same-rf-pairs", type=int, default=2)
    p.add_argument("--drop", type=float, default=0.10)
    p.add_argument("--spurious", type=float, default=0.02)
    p.add_argument("--library", help="emitter library CSV; also scores identification")
    p.add_argument("--unknown-frac", type=float, default=0.25)
    a = p.parse_args()

    found = collections.Counter()
    total = collections.Counter()
    ident = collections.Counter()
    n_out = n_true = 0
    with tempfile.TemporaryDirectory() as work:
        for seed in range(a.runs):
            f, t, i, o, n = score_one(a.cli, seed, a.emitters, a.same_rf_pairs, a.drop,
                                      a.spurious, work, a.library, a.unknown_frac)
            found.update(f)
            total.update(t)
            ident.update(i)
            n_out += o
            n_true += n

    print(f"{a.runs} scenarios, {a.emitters} emitters each, drop {a.drop:.0%}, spurious {a.spurious:.0%}")
    for kind in ("constant", "jitter", "stagger"):
        print(f"  {kind:<9} {found[kind]:4d} / {total[kind]:<4d} found")
    print(f"  emitters reported: {n_out} (true: {n_true})")
    if a.library:
        print("identification (found emitters only):")
        for k in ("known, right name", "known, called UNKNOWN", "known, wrong name",
                  "unknown, called UNKNOWN", "unknown, given a name"):
            print(f"  {k:<24} {ident[k]:4d}")


if __name__ == "__main__":
    main()
