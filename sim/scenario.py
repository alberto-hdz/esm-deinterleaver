"""Generate a synthetic ESM/RWR scenario as a list of PDWs.

I simulate a handful of made-up emitters, build each one's pulse train,
add receiver measurement noise, drop some pulses, add some spurious ones,
and write everything out sorted by time of arrival.

Outputs (in --out):
  pdw.csv   toa_us,freq_mhz,pw_us,amp_db,aoa_deg,truth_id   (truth_id -1 = spurious)
  truth.csv one row per emitter with its true parameters
"""
import argparse
import csv
from pathlib import Path

import numpy as np

# Receiver measurement noise (1 sigma)
TOA_SIGMA_US = 0.05
RF_SIGMA_MHZ = 2.0
PW_SIGMA_US = 0.05
AOA_SIGMA_DEG = 3.0
AMP_SIGMA_DB = 1.0


def pulse_train(pri, kind, T, rng, jitter=0.0, stagger=None):
    """Return pulse times of arrival (us) for one emitter over a T us window."""
    n = int(T / (min(stagger) if stagger else pri)) + 1
    if kind == "constant":
        iv = np.full(n, pri)
    elif kind == "jitter":
        iv = pri * (1 + rng.uniform(-jitter, jitter, n))
    elif kind == "stagger":
        iv = np.resize(stagger, n)
    else:
        raise ValueError(f"unknown PRI type {kind}")
    # random start phase so the trains don't all line up at t = 0
    toa = rng.uniform(0, iv[0]) + np.concatenate([[0], np.cumsum(iv[:-1])])
    return toa[toa < T]


def make_emitter(eid, rng):
    pri = float(np.exp(rng.uniform(np.log(100), np.log(3000))))
    pw = float(np.exp(rng.uniform(np.log(0.2), np.log(50))))
    pw = min(pw, 0.1 * pri)  # keep duty cycle at or below 10%
    kind = rng.choice(["constant", "jitter", "stagger"], p=[0.5, 0.3, 0.2])
    jitter, stagger = 0.0, None
    if kind == "jitter":
        jitter = float(rng.uniform(0.01, 0.10))
    elif kind == "stagger":
        levels = int(rng.integers(2, 5))
        stagger = [round(pri * (1 + rng.uniform(-0.3, 0.3)), 1) for _ in range(levels)]
        pri = float(np.mean(stagger))
    return {
        "id": eid,
        "freq_mhz": float(rng.uniform(2000, 18000)),
        "pw_us": pw,
        "pri_type": str(kind),
        "pri_us": pri,
        "jitter": jitter,
        "stagger": stagger,
        "aoa_deg": float(rng.uniform(0, 360)),
        "amp_db": float(rng.uniform(-60, -30)),
    }


def build_scenario(n_emitters, drop, spurious, duration_ms, same_rf_pairs, seed):
    rng = np.random.default_rng(seed)
    T = duration_ms * 1000.0
    emitters = [make_emitter(i, rng) for i in range(n_emitters)]

    # Hard case: copy RF, AOA and PW from emitter 2k onto emitter 2k+1 so the
    # clusterer can't split them. Only PRI analysis can separate these pairs.
    for k in range(min(same_rf_pairs, n_emitters // 2)):
        a, b = emitters[2 * k], emitters[2 * k + 1]
        b["freq_mhz"] = a["freq_mhz"] + rng.uniform(-1, 1)
        b["aoa_deg"] = (a["aoa_deg"] + rng.uniform(-2, 2)) % 360
        b["pw_us"] = min(a["pw_us"], 0.1 * b["pri_us"])

    cols = {k: [] for k in ["toa_us", "freq_mhz", "pw_us", "amp_db", "aoa_deg", "truth_id"]}
    for e in emitters:
        toa = pulse_train(e["pri_us"], e["pri_type"], T, rng, e["jitter"], e["stagger"])
        keep = rng.random(toa.size) >= drop
        toa = toa[keep]
        n = toa.size
        e["n_pulses"] = n
        cols["toa_us"].append(toa + rng.normal(0, TOA_SIGMA_US, n))
        cols["freq_mhz"].append(e["freq_mhz"] + rng.normal(0, RF_SIGMA_MHZ, n))
        cols["pw_us"].append(np.maximum(e["pw_us"] + rng.normal(0, PW_SIGMA_US, n), 0.05))
        cols["amp_db"].append(e["amp_db"] + rng.normal(0, AMP_SIGMA_DB, n))
        cols["aoa_deg"].append((e["aoa_deg"] + rng.normal(0, AOA_SIGMA_DEG, n)) % 360)
        cols["truth_id"].append(np.full(n, e["id"]))

    # Spurious pulses: a fraction of the real pulse count with random parameters
    n_real = sum(len(c) for c in cols["toa_us"])
    n_sp = int(round(spurious * n_real))
    cols["toa_us"].append(rng.uniform(0, T, n_sp))
    cols["freq_mhz"].append(rng.uniform(2000, 18000, n_sp))
    cols["pw_us"].append(np.exp(rng.uniform(np.log(0.2), np.log(50), n_sp)))
    cols["amp_db"].append(rng.uniform(-70, -30, n_sp))
    cols["aoa_deg"].append(rng.uniform(0, 360, n_sp))
    cols["truth_id"].append(np.full(n_sp, -1))

    pdw = {k: np.concatenate(v) for k, v in cols.items()}
    order = np.argsort(pdw["toa_us"])
    pdw = {k: v[order] for k, v in pdw.items()}
    return emitters, pdw


def write_outputs(out, emitters, pdw):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    with open(out / "pdw.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["toa_us", "freq_mhz", "pw_us", "amp_db", "aoa_deg", "truth_id"])
        for i in range(pdw["toa_us"].size):
            w.writerow([f"{pdw['toa_us'][i]:.4f}", f"{pdw['freq_mhz'][i]:.3f}",
                        f"{pdw['pw_us'][i]:.4f}", f"{pdw['amp_db'][i]:.2f}",
                        f"{pdw['aoa_deg'][i]:.2f}", int(pdw["truth_id"][i])])
    with open(out / "truth.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["id", "freq_mhz", "pw_us", "pri_type", "pri_us", "jitter",
                    "stagger_us", "aoa_deg", "amp_db", "n_pulses"])
        for e in emitters:
            stag = ";".join(str(s) for s in e["stagger"]) if e["stagger"] else ""
            w.writerow([e["id"], f"{e['freq_mhz']:.3f}", f"{e['pw_us']:.4f}", e["pri_type"],
                        f"{e['pri_us']:.3f}", f"{e['jitter']:.4f}", stag,
                        f"{e['aoa_deg']:.2f}", f"{e['amp_db']:.2f}", e["n_pulses"]])


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--emitters", type=int, default=6)
    p.add_argument("--drop", type=float, default=0.10, help="fraction of pulses lost")
    p.add_argument("--spurious", type=float, default=0.02, help="spurious pulses as a fraction of real ones")
    p.add_argument("--duration-ms", type=float, default=200.0)
    p.add_argument("--same-rf-pairs", type=int, default=1, help="emitter pairs sharing RF/AOA/PW")
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--out", required=True)
    a = p.parse_args()

    emitters, pdw = build_scenario(a.emitters, a.drop, a.spurious, a.duration_ms,
                                   a.same_rf_pairs, a.seed)
    write_outputs(a.out, emitters, pdw)

    print(f"{pdw['toa_us'].size} PDWs ({(pdw['truth_id'] < 0).sum()} spurious) -> {a.out}")
    for e in emitters:
        extra = f" jitter={e['jitter']:.1%}" if e["pri_type"] == "jitter" else ""
        extra += f" stagger={e['stagger']}" if e["stagger"] else ""
        print(f"  E{e['id']}: RF {e['freq_mhz']:8.1f} MHz  PW {e['pw_us']:6.2f} us  "
              f"PRI {e['pri_us']:7.1f} us ({e['pri_type']}{extra})  AOA {e['aoa_deg']:5.1f}  "
              f"pulses {e['n_pulses']}")


if __name__ == "__main__":
    main()
