"""Make the result plots for the README.

Three Monte Carlo sweeps, all built on the same scoring as score.py and score_iq.py:
  docs/found_vs_loss.png      emitters found vs pulse loss, by PRI type
  docs/found_vs_emitters.png  emitters found vs how many emitters are on the air
  docs/iq_vs_snr.png          pulse detection and TOA error vs SNR (raw I/Q path)

usage: python sim/plots.py --cli build/esm_cli --out docs
Takes a few minutes; --runs sets the scenarios per point.
"""
import argparse
import collections
import sys
import tempfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import score  # noqa: E402
import score_iq  # noqa: E402

# Light theme. The three PRI types always get the same color in every plot.
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_2 = "#52514e"
GRID = "#e4e3df"
COLORS = {"constant": "#2a78d6", "jitter": "#eb6834", "stagger": "#1baf7a"}
MARKERS = {"constant": "o", "jitter": "s", "stagger": "^"}


def style(ax, title, xlabel, ylabel):
    ax.set_facecolor(SURFACE)
    ax.set_title(title, loc="left", color=INK, fontsize=11, pad=10)
    ax.set_xlabel(xlabel, color=INK_2, fontsize=9)
    ax.set_ylabel(ylabel, color=INK_2, fontsize=9)
    ax.grid(axis="y", color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)
    ax.tick_params(colors=INK_2, labelsize=8, length=0)


def by_type_plot(path, title, xlabel, xs, pct):
    fig, ax = plt.subplots(figsize=(6.4, 3.6), dpi=150, facecolor=SURFACE)
    style(ax, title, xlabel, "emitters found (%)")
    for kind, ys in pct.items():
        ax.plot(xs, ys, color=COLORS[kind], linewidth=1.5, marker=MARKERS[kind], markersize=5,
                markeredgecolor=SURFACE, markeredgewidth=1.2, label=kind)

    # Direct labels at the end of each line, in text ink (the marker carries the color).
    # Lines often end at the same value, so I push the labels apart to keep them readable.
    ends = sorted(((ys[-1], kind) for kind, ys in pct.items()), reverse=True)
    label_y = []
    for y, _ in ends:
        label_y.append(min(y, label_y[-1] - 5) if label_y else y)
    for (y, kind), ly in zip(ends, label_y):
        ax.annotate(kind, (xs[-1], ly), xytext=(8, 0), textcoords="offset points",
                    va="center", fontsize=8, color=INK_2)
    ax.set_xticks(xs)
    ax.set_ylim(0, 105)
    ax.set_xlim(xs[0], xs[-1] + (xs[-1] - xs[0]) * 0.15)
    ax.legend(frameon=False, fontsize=8, labelcolor=INK_2, loc="lower left")
    fig.tight_layout()
    fig.savefig(path, facecolor=SURFACE)
    plt.close(fig)


def sweep(cli, runs, values, make_args):
    pct = {k: [] for k in COLORS}
    with tempfile.TemporaryDirectory() as work:
        for v in values:
            found = collections.Counter()
            total = collections.Counter()
            emitters, drop = make_args(v)
            for seed in range(runs):
                f, t, _, _, _ = score.score_one(cli, seed, emitters, 2, drop, 0.02, work, None, 0.25)
                found.update(f)
                total.update(t)
            for k in COLORS:
                pct[k].append(100.0 * found[k] / max(total[k], 1))
            print(f"  {v}: " + "  ".join(f"{k} {found[k]}/{total[k]}" for k in COLORS))
    return pct


def snr_sweep(cli, runs, snrs):
    pd, toa = [], []
    with tempfile.TemporaryDirectory() as work:
        for snr in snrs:
            clean = det = 0
            errs = []
            for seed in range(runs):
                s, _, _ = score_iq.score_one(cli, seed, 1, work,
                                             ["--snr-db", str(snr), "--duration-ms", "20"])
                clean += s["clean"]
                det += s["clean_detected"]
                errs.append(s["toa_err_ns"])
            errs = np.concatenate(errs)
            pd.append(100.0 * det / clean)
            toa.append(float(np.sqrt(np.mean(errs ** 2))) if errs.size else np.nan)
            print(f"  {snr:4.1f} dB: detected {det}/{clean}, TOA RMS {toa[-1]:.1f} ns")
    return pd, toa


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--cli", default="build/esm_cli")
    p.add_argument("--out", default="docs")
    p.add_argument("--runs", type=int, default=20)
    a = p.parse_args()
    out = Path(a.out)
    out.mkdir(exist_ok=True)

    print("pulse loss sweep")
    loss = [0, 10, 20, 30, 40]
    pct = sweep(a.cli, a.runs, loss, lambda v: (8, v / 100))
    by_type_plot(out / "found_vs_loss.png", "Emitters found vs pulse loss (8 emitters)",
                 "pulses lost (%)", loss, pct)

    print("emitter count sweep")
    counts = [4, 8, 12, 16]
    pct = sweep(a.cli, a.runs, counts, lambda v: (v, 0.10))
    by_type_plot(out / "found_vs_emitters.png", "Emitters found vs emitters on the air (10% loss)",
                 "emitters in the scenario", counts, pct)

    print("I/Q SNR sweep")
    snrs = [0, 2, 4, 6, 8, 10, 12, 16, 20]
    pd, toa = snr_sweep(a.cli, max(a.runs // 4, 3), snrs)
    # two separate panels, each with its own y axis (never two scales on one plot)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.2, 3.2), dpi=150, facecolor=SURFACE)
    style(ax1, "Pulses detected", "SNR per sample (dB)", "detected (%)")
    ax1.plot(snrs, pd, color=COLORS["constant"], linewidth=1.5, marker="o", markersize=5,
             markeredgecolor=SURFACE, markeredgewidth=1.2)
    ax1.set_ylim(0, 105)
    style(ax2, "TOA error (RMS)", "SNR per sample (dB)", "ns  (sample period = 25 ns)")
    keep = [i for i, v in enumerate(pd) if v >= 50]  # TOA error only means something once pulses are found
    ax2.plot([snrs[i] for i in keep], [toa[i] for i in keep], color=COLORS["constant"],
             linewidth=1.5, marker="o", markersize=5, markeredgecolor=SURFACE, markeredgewidth=1.2)
    ax2.set_ylim(bottom=0)
    for ax in (ax1, ax2):
        ax.set_xlim(-1, 21)
        ax.set_xticks([0, 5, 10, 15, 20])
    fig.tight_layout()
    fig.savefig(out / "iq_vs_snr.png", facecolor=SURFACE)
    plt.close(fig)
    print(f"wrote plots to {out}/")


if __name__ == "__main__":
    main()
