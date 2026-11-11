import csv
import subprocess
import sys
from pathlib import Path

import numpy as np

GEN = Path(__file__).resolve().parents[1] / "iq_gen.py"


def run_gen(tmp_path, *args):
    subprocess.run([sys.executable, str(GEN), "--out", str(tmp_path), *args],
                   check=True, capture_output=True)
    x = np.fromfile(tmp_path / "iq.c64", dtype=np.complex64)
    pulses = list(csv.DictReader(open(tmp_path / "truth_pulses.csv")))
    return x, pulses


def test_size_and_noise_power(tmp_path):
    x, _ = run_gen(tmp_path, "--duration-ms", "5", "--fs-mhz", "40", "--seed", "2")
    assert x.size == 200_000
    # median-based like the C++ detector, so the pulses don't skew it
    noise = np.median(np.abs(x) ** 2) / np.log(2)
    assert 0.9 < noise < 1.2


def test_pulse_lands_where_truth_says(tmp_path):
    x, pulses = run_gen(tmp_path, "--duration-ms", "5", "--seed", "3", "--emitters", "1")
    p = pulses[len(pulses) // 2]
    k = int((float(p["toa_us"]) + float(p["pw_us"]) / 2) * 40)  # middle of the pulse
    in_pulse = np.mean(np.abs(x[k - 5:k + 5]) ** 2)
    assert in_pulse > 5.0  # SNR is at least 8 dB, so well above the unit noise


def test_capture_file_has_rates(tmp_path):
    run_gen(tmp_path, "--duration-ms", "1", "--fs-mhz", "20", "--fc-mhz", "3000")
    assert (tmp_path / "capture.txt").read_text().split() == ["--fs-mhz", "20", "--fc-mhz", "3000"]
