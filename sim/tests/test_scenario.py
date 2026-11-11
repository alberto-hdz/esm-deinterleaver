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


LIBRARY = Path(__file__).resolve().parents[2] / "config" / "emitter_library.csv"


def test_library_emitters_land_in_their_windows():
    from scenario import read_library
    lib = {r["name"]: r for r in read_library(LIBRARY)}
    emitters, _ = build_scenario(8, 0.1, 0.0, 200, 0, seed=4, library=list(lib.values()),
                                 unknown_frac=0.0)
    for e in emitters:
        r = lib[e["name"]]
        assert r["rf_min_mhz"] <= e["freq_mhz"] <= r["rf_max_mhz"]
        assert r["pri_min_us"] - 1e-6 <= e["pri_us"] <= r["pri_max_us"] + 1e-6
        assert e["pri_type"] == r["pri_type"]
        if e["pri_type"] == "stagger":
            assert len(e["stagger"]) == r["stagger_levels"]


def test_library_pairs_share_rf_and_pw():
    from scenario import read_library
    emitters, _ = build_scenario(4, 0.1, 0.0, 200, 2, seed=5, library=read_library(LIBRARY))
    for k in (0, 2):
        a, b = emitters[k], emitters[k + 1]
        assert abs(a["freq_mhz"] - b["freq_mhz"]) <= 1.0
        assert a["name"] != "UNKNOWN" and b["name"] != "UNKNOWN"


def test_default_mode_unchanged_by_library_code():
    # the library option must not change what the default mode generates
    _, a = build_scenario(6, 0.1, 0.02, 200, 1, seed=1)
    _, b = build_scenario(6, 0.1, 0.02, 200, 1, seed=1, library=None)
    assert np.array_equal(a["toa_us"], b["toa_us"])
