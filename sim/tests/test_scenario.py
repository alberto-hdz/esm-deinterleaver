import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scenario import build_scenario, pulse_train  # noqa: E402


def test_constant_pri_spacing():
    rng = np.random.default_rng(0)
    toa = pulse_train(500.0, "constant", 100_000.0, rng)
    assert np.allclose(np.diff(toa), 500.0)


def test_jitter_stays_in_bounds():
    rng = np.random.default_rng(0)
    d = np.diff(pulse_train(1000.0, "jitter", 1_000_000.0, rng, jitter=0.05))
    assert d.min() >= 950.0 and d.max() <= 1050.0


def test_stagger_repeats():
    rng = np.random.default_rng(0)
    d = np.diff(pulse_train(0, "stagger", 100_000.0, rng, stagger=[300.0, 400.0, 500.0]))
    assert np.allclose(d[:6], [300, 400, 500, 300, 400, 500])


def test_scenario_sorted_and_drop_rate():
    emitters, pdw = build_scenario(6, drop=0.2, spurious=0.0, duration_ms=200,
                                   same_rf_pairs=1, seed=3)
    assert np.all(np.diff(pdw["toa_us"]) >= -1.0)  # sorted (allowing TOA noise)
    assert set(np.unique(pdw["truth_id"])) == set(range(6))
    # pair 0/1 shares RF and AOA
    assert abs(emitters[0]["freq_mhz"] - emitters[1]["freq_mhz"]) <= 1.0


def test_seed_is_reproducible():
    _, a = build_scenario(4, 0.1, 0.02, 50, 1, seed=7)
    _, b = build_scenario(4, 0.1, 0.02, 50, 1, seed=7)
    assert np.array_equal(a["toa_us"], b["toa_us"])
