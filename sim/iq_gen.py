"""Generate raw complex I/Q samples from a receiver watching a few emitters.

This is one step earlier than scenario.py: instead of a ready-made list of PDWs,
it writes what a single-channel receiver would digitize, so the C++ pulse
detector has to find the pulses and measure them itself.

Every emitter sits at a frequency offset inside the receiver's band and sends
unmodulated pulses with 25 ns edges. Noise is complex Gaussian with unit power,
so an emitter's SNR in dB is just 10*log10(amplitude^2).

Outputs (in --out):
  iq.c64             raw complex64 samples (interleaved float32 I, Q)
  truth_pulses.csv   toa_us,freq_mhz,pw_us,truth_id   one row per transmitted pulse
  truth.csv          one row per emitter
  capture.txt        fs and fc, for the esm_cli command line
"""
import argparse
import csv
from pathlib import Path

import numpy as np

from scenario import pulse_train

EDGE_S = 25e-9


def make_emitters(n, fs_mhz, rng, same_rf_pairs):
    # Frequency slots 12 MHz apart, more than the clusterer's 10 MHz tolerance.
    # With more emitters than slots, some share a slot and only PW or PRI can
    # separate them, same idea as the forced pairs in scenario.py.
    half_band = fs_mhz / 2 - 3
    slots = np.arange(-half_band, half_band + 0.1, 12.0)
    rng.shuffle(slots)
    emitters = []
    for i in range(n):
        pri = float(np.exp(rng.uniform(np.log(100), np.log(1000))))
        kind = str(rng.choice(["constant", "jitter", "stagger"], p=[0.5, 0.3, 0.2]))
        jitter, stagger = 0.0, None
        if kind == "jitter":
            jitter = float(rng.uniform(0.01, 0.10))
        elif kind == "stagger":
            iv = pri * (1 + rng.uniform(-0.3, 0.3, int(rng.integers(2, 5))))
            stagger = [round(float(v), 1) for v in iv]
            pri = float(np.mean(stagger))
        emitters.append({
            "id": i,
            "offset_mhz": float(slots[i % len(slots)] + rng.uniform(-1, 1)),
            # keep the duty cycle at or below 3% so pulses from different emitters
            # rarely overlap; overlapping pulses merge in a single-channel receiver
            "pw_us": float(min(rng.uniform(0.5, 10.0), 0.03 * min(stagger or [pri]))),
            "pri_type": kind,
            "pri_us": pri,
            "jitter": jitter,
            "stagger": stagger,
            "snr_db": float(rng.uniform(8, 25)),
        })
    for k in range(min(same_rf_pairs, n // 2)):
        a, b = emitters[2 * k], emitters[2 * k + 1]
        b["offset_mhz"] = a["offset_mhz"] + rng.uniform(-0.5, 0.5)
        b["pw_us"] = min(a["pw_us"], 0.03 * min(b["stagger"] or [b["pri_us"]]))
    return emitters


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--emitters", type=int, default=4)
    p.add_argument("--duration-ms", type=float, default=50.0)
    p.add_argument("--fs-mhz", type=float, default=40.0)
    p.add_argument("--fc-mhz", type=float, default=9000.0)
    p.add_argument("--same-rf-pairs", type=int, default=1)
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--out", required=True)
    a = p.parse_args()

    rng = np.random.default_rng(a.seed)
    fs = a.fs_mhz * 1e6
    n = int(a.duration_ms * 1e-3 * fs)
    T_us = a.duration_ms * 1000.0
    emitters = make_emitters(a.emitters, a.fs_mhz, rng, a.same_rf_pairs)

    # complex Gaussian noise, unit power per sample
    x = (rng.standard_normal(n) + 1j * rng.standard_normal(n)) * np.sqrt(0.5)

    truth_rows = []
    for e in emitters:
        toa = pulse_train(e["pri_us"], e["pri_type"], T_us, rng, e["jitter"], e["stagger"])
        e["n_pulses"] = int(toa.size)
        amp = np.sqrt(10 ** (e["snr_db"] / 10))
        f = e["offset_mhz"] * 1e6
        for t0_us in toa:
            t0 = t0_us * 1e-6
            t1 = t0 + e["pw_us"] * 1e-6
            k0 = max(int((t0 - EDGE_S) * fs), 0)
            k1 = min(int((t1 + EDGE_S) * fs) + 2, n)
            t = np.arange(k0, k1) / fs
            # linear 25 ns edges, evaluated at the real sample times so a pulse can
            # start between samples; TOA and PW are measured at the 50% points
            env = np.clip((t - t0) / EDGE_S + 0.5, 0, 1) * np.clip((t1 - t) / EDGE_S + 0.5, 0, 1)
            x[k0:k1] += amp * env * np.exp(2j * np.pi * f * t)
            truth_rows.append((t0_us, a.fc_mhz + e["offset_mhz"], e["pw_us"], e["id"]))

    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    x.astype(np.complex64).tofile(out / "iq.c64")
    truth_rows.sort()
    with open(out / "truth_pulses.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["toa_us", "freq_mhz", "pw_us", "truth_id"])
        for r in truth_rows:
            w.writerow([f"{r[0]:.4f}", f"{r[1]:.3f}", f"{r[2]:.4f}", r[3]])
    with open(out / "truth.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["id", "freq_mhz", "pw_us", "pri_type", "pri_us", "jitter", "stagger_us",
                    "snr_db", "n_pulses"])
        for e in emitters:
            stag = ";".join(str(s) for s in e["stagger"]) if e["stagger"] else ""
            w.writerow([e["id"], f"{a.fc_mhz + e['offset_mhz']:.3f}", f"{e['pw_us']:.4f}",
                        e["pri_type"], f"{e['pri_us']:.3f}", f"{e['jitter']:.4f}", stag,
                        f"{e['snr_db']:.1f}", e["n_pulses"]])
    (out / "capture.txt").write_text(f"--fs-mhz {a.fs_mhz:g} --fc-mhz {a.fc_mhz:g}\n")

    print(f"{n} samples ({a.duration_ms:g} ms at {a.fs_mhz:g} MSPS, {n * 8 / 1e6:.0f} MB), "
          f"{len(truth_rows)} pulses -> {out}")
    for e in emitters:
        extra = f" jitter={e['jitter']:.1%}" if e["pri_type"] == "jitter" else ""
        extra += f" stagger={e['stagger']}" if e["stagger"] else ""
        print(f"  E{e['id']}: RF {a.fc_mhz + e['offset_mhz']:8.1f} MHz  PW {e['pw_us']:5.2f} us  "
              f"PRI {e['pri_us']:6.1f} us ({e['pri_type']}{extra})  SNR {e['snr_db']:4.1f} dB  "
              f"pulses {e['n_pulses']}")


if __name__ == "__main__":
    main()
