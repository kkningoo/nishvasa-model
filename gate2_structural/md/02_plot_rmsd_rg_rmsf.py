#!/usr/bin/env python3
"""
Plot RMSD, radius of gyration, and per-residue RMSF from gmx .xvg output.

RMSF is split into chain A / chain B (two stacked panels + an overlay) rather
than one continuous plot -- for a homodimer, gmx rmsf's per-residue output is
indexed continuously across both chains (e.g. rows 1-443 = chain A, rows
444-886 = chain B, each chain internally renumbered from 1), so a single
unbroken line is confusing to read and hides a very real question: do the two
chains behave the same?

Usage:
    python 02_plot_rmsd_rg_rmsf.py <rmsd.xvg> <gyrate.xvg> <rmsf.xvg> <out_prefix> <title> [chainA_length]

chainA_length defaults to 443 (the GS dimer construct used in this project --
confirmed by checking where gmx rmsf's printed residue-number column restarts
at 1, e.g. `grep -v '^[#@]' rmsf.xvg | awk 'NR==443||NR==444'`). Pass a
different value for other systems.
"""
import sys
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def load_xvg(path):
    t, y = [], []
    with open(path) as f:
        for line in f:
            if line.startswith(("#", "@")):
                continue
            parts = line.split()
            if len(parts) >= 2:
                t.append(float(parts[0]))
                y.append(float(parts[1]))
    return t, y


def safe_save(fig, path):
    """Write to a temp name then rename over the target -- on some
    WSL/drvfs-mounted Windows drives, overwriting a PNG in place with
    matplotlib's savefig can fail with OSError (EINVAL); this sidesteps it."""
    tmp = path + ".tmp.png"
    fig.savefig(tmp, dpi=120)
    plt.close(fig)
    os.replace(tmp, path)


def plot_rmsd(xvg_path, out_png, title):
    t, y = load_xvg(xvg_path)
    fig = plt.figure(figsize=(7, 4))
    plt.plot(t, y)
    plt.xlabel("Time (ns)")
    plt.ylabel("RMSD (nm)")
    plt.title(f"{title} -- backbone RMSD")
    plt.tight_layout()
    safe_save(fig, out_png)
    print("wrote", out_png)


def plot_gyrate(xvg_path, out_png, title):
    t, y = load_xvg(xvg_path)
    fig = plt.figure(figsize=(7, 4))
    plt.plot(t, y)
    plt.xlabel("Time (ps)")
    plt.ylabel("Rg (nm)")
    plt.title(f"{title} -- radius of gyration")
    plt.tight_layout()
    safe_save(fig, out_png)
    print("wrote", out_png)


def plot_rmsf_by_chain(xvg_path, out_png, title, chainA_length=443):
    _, vals = load_xvg(xvg_path)
    chainA = vals[:chainA_length]
    chainB = vals[chainA_length:]
    xA = list(range(1, len(chainA) + 1))
    xB = list(range(1, len(chainB) + 1))

    fig, axes = plt.subplots(2, 1, figsize=(9, 7), sharex=True)
    axes[0].plot(xA, chainA, color="tab:blue")
    axes[0].set_title(f"{title} -- Chain A")
    axes[0].set_ylabel("RMSF (nm)")
    axes[1].plot(xB, chainB, color="tab:orange")
    axes[1].set_title(f"{title} -- Chain B")
    axes[1].set_ylabel("RMSF (nm)")
    axes[1].set_xlabel("Residue (local chain numbering)")
    fig.tight_layout()
    safe_save(fig, out_png)
    print("wrote", out_png, "chainA n=", len(chainA), "chainB n=", len(chainB))

    fig2 = plt.figure(figsize=(9, 4))
    plt.plot(xA, chainA, label="Chain A", color="tab:blue", lw=0.9)
    plt.plot(xB, chainB, label="Chain B", color="tab:orange", lw=0.9)
    plt.xlabel("Residue (local chain numbering)")
    plt.ylabel("RMSF (nm)")
    plt.title(f"{title} -- Chain A vs Chain B overlay")
    plt.legend()
    plt.tight_layout()
    overlay_png = out_png.replace(".png", "_overlay.png")
    safe_save(fig2, overlay_png)
    print("wrote", overlay_png)


if __name__ == "__main__":
    if len(sys.argv) < 6:
        print(__doc__)
        sys.exit(1)
    rmsd_xvg, gyrate_xvg, rmsf_xvg, out_prefix, title = sys.argv[1:6]
    chainA_len = int(sys.argv[6]) if len(sys.argv) > 6 else 443

    plot_rmsd(rmsd_xvg, out_prefix + "_rmsd.png", title)
    plot_gyrate(gyrate_xvg, out_prefix + "_gyrate.png", title)
    plot_rmsf_by_chain(rmsf_xvg, out_prefix + "_rmsf.png", title, chainA_len)
