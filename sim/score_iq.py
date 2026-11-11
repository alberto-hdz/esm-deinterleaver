"""Score the I/Q path: pulse detection accuracy, then emitters found end to end.

For each seed I generate a capture with iq_gen.py and run esm_cli --iq on it.
Each transmitted pulse is matched to the detected pulse closest in time (within
0.5 us). From the matches I get probability of detection, false alarms, and the
TOA / PW / RF measurement errors. Pulses that overlap another emitter's pulse in
time are counted separately: a single-channel receiver can't pull those apart, so
they come out as one merged pulse with the wrong width and frequency. Then every
true emitter is checked the same way score.py does it, using the matched pulses.

usage: python sim/score_iq.py --cli build/esm_cli --runs 20
"""
import argparse
import collections
import csv
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

GEN = Path(__file__).resolve().parent / "iq_gen.py"
MATCH_US = 0.5


def load(path, *cols):
    rows = list(csv.DictReader(open(path)))
    return [np.array([float(r[c]) for r in rows]) for c in cols]


def score_one(cli, seed, emitters, workdir):
    out = Path(workdir) / f"iq{seed}"
    subprocess.run([sys.executable, str(GEN), "--emitters", str(emitters), "--seed", str(seed),
                    "--out", str(out)], check=True, capture_output=True)
    cap = (out / "capture.txt").read_text().split()
    subprocess.run([cli, "--iq", str(out / "iq.c64"), *cap, "--pdw-out", str(out / "det.csv"),
                    "--labels", str(out / "labels.csv")], check=True, capture_output=True)

    t_toa, t_rf, t_pw, t_id = load(out / "truth_pulses.csv", "toa_us", "freq_mhz", "pw_us", "truth_id")
    d_toa, d_rf, d_pw = load(out / "det.csv", "toa_us", "freq_mhz", "pw_us")
    (lab,) = load(out / "labels.csv", "emitter")

    # nearest detection for every true pulse (both lists are sorted by time)
    j = np.clip(np.searchsorted(d_toa, t_toa), 1, max(len(d_toa) - 1, 1))
    left_closer = np.abs(t_toa - d_toa[j - 1]) < np.abs(t_toa - d_toa[j])
    j = np.where(left_closer, j - 1, j)
    hit = np.abs(t_toa - d_toa[j]) <= MATCH_US
    used = np.zeros(len(d_toa), bool)
    used[j[hit]] = True

    # a pulse "overlaps" if any other pulse is on the air within 0.3 us of it
    t_end = t_toa + t_pw
    order = np.argsort(t_toa)
    overlap = np.zeros(len(t_toa), bool)
    for a_, b_ in zip(order[:-1], order[1:]):
        if t_toa[b_] < t_end[a_] + 0.3:
            overlap[a_] = overlap[b_] = True
    clean = hit & ~overlap

    stats = {
        "pulses": len(t_toa),
        "overlapped": int(overlap.sum()),
        "clean": int((~overlap).sum()),
        "clean_detected": int(clean.sum()),
        # detections that matched no pulse: either a leftover piece of an overlap, or
        # a real false alarm in noise-only time
        "extra_in_pulse": int(sum(np.any((t_toa - MATCH_US <= x) & (x <= t_end + MATCH_US))
                                  for x in d_toa[~used])),
        "false": int((~used).sum()),
        "toa_err_ns": (d_toa[j[clean]] - t_toa[clean]) * 1000,
        "pw_err_ns": (d_pw[j[clean]] - t_pw[clean]) * 1000,
        "rf_err_khz": (d_rf[j[clean]] - t_rf[clean]) * 1000,
    }

    # end to end: give each detection the truth id of the pulse it matched, then
    # count an emitter as found if one output emitter has >= 80% of its detected
    # pulses at >= 90% purity
    det_id = np.full(len(d_toa), -1)
    det_id[j[hit]] = t_id[hit].astype(int)
    truth = list(csv.DictReader(open(out / "truth.csv")))
    found = collections.Counter()
    total = collections.Counter()
    for t in truth:
        mine = det_id == int(t["id"])
        total[t["pri_type"]] += 1
        if mine.sum() == 0:
            continue
        labels, counts = np.unique(lab[mine], return_counts=True)
        for L, cnt in zip(labels, counts):
            if L >= 0 and cnt / mine.sum() >= 0.8 and cnt / (lab == L).sum() >= 0.9:
                found[t["pri_type"]] += 1
                break
    return stats, found, total


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--cli", default="build/esm_cli")
    p.add_argument("--runs", type=int, default=20)
    p.add_argument("--emitters", type=int, default=4)
    a = p.parse_args()

    pulses = overlapped = clean = clean_det = false = extra = 0
    toa, pw, rf = [], [], []
    found = collections.Counter()
    total = collections.Counter()
    with tempfile.TemporaryDirectory() as work:
        for seed in range(a.runs):
            s, f, t = score_one(a.cli, seed, a.emitters, work)
            pulses += s["pulses"]
            overlapped += s["overlapped"]
            clean += s["clean"]
            clean_det += s["clean_detected"]
            false += s["false"]
            extra += s["extra_in_pulse"]
            toa.append(s["toa_err_ns"])
            pw.append(s["pw_err_ns"])
            rf.append(s["rf_err_khz"])
            found.update(f)
            total.update(t)
    toa, pw, rf = np.concatenate(toa), np.concatenate(pw), np.concatenate(rf)
    rms = lambda v: float(np.sqrt(np.mean(v ** 2)))

    print(f"{a.runs} captures, {a.emitters} emitters each, SNR 8-25 dB")
    print(f"  pulses transmitted        {pulses}")
    print(f"  overlapping another pulse {overlapped} ({overlapped / pulses:.1%}), not scored below")
    print(f"  clean pulses detected     {clean_det} / {clean} ({clean_det / clean:.2%})")
    print(f"  extra detections during real pulses (overlap pieces) {extra}")
    print(f"  false alarms in noise-only time                      {false - extra}")
    print("measurement error on clean pulses:")
    print(f"  TOA  RMS {rms(toa):6.1f} ns    (sample period is 25 ns)")
    print(f"  PW   RMS {rms(pw):6.1f} ns")
    print(f"  RF   RMS {rms(rf):6.1f} kHz")
    print("emitters found end to end:")
    for kind in ("constant", "jitter", "stagger"):
        print(f"  {kind:<9} {found[kind]:4d} / {total[kind]:<4d}")


if __name__ == "__main__":
    main()
