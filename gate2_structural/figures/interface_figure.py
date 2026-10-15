"""Interface energetics, two panels.

(a) energy density per interface, with buried area shown and the two runs
    that did not converge drawn rather than omitted
(b) the stitch-count comparison that this result reversed

    python3 interface_figure.py

Writes interface_energetics.png and, with --split, two separate files.
"""

import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch

# ---- inputs --------------------------------------------------------------
# dG per 1000 A^2 is reported positive here; the axis is labelled -REU so
# taller means stronger. buried_area is the total across both sides.
#          label            dG/1000A2  buried  group   converged
BARS = [
    ('A–H\ncross-ring',        21.20,   4251,  'gs',   True),
    ('A–F\nintra-ring',        12.85,   3460,  'gs',   True),
    ('A–B\nintra-ring',        10.28,   3477,  'gs',   True),
    ('C–D\nintra-ring',          None,  3316,  'gs',   False),
    ('A–C',                    19.96,   2029,  'ppk2', True),
    ('A–B',                    18.34,   1795,  'ppk2', True),
    ('C–D',                      None,  1846,  'ppk2', False),
]

# genipin span 5-16 A, cross-subunit pairs only
STITCH = {'within-ring': 26, 'cross-ring': 6}

COL = {'gs_strong': '#534AB7', 'gs_weak': '#888780',
       'ppk2': '#1D9E75', 'open': '#B0ADA6'}


def panel_a(ax):
    x = np.arange(len(BARS))
    labels, heights, areas, cols, hatches, edges = [], [], [], [], [], []

    for lab, dg, area, grp, conv in BARS:
        labels.append(lab)
        areas.append(area)
        if not conv:
            # draw at the panel's mean so the bar is visible, but unfilled
            heights.append(14.0)
            cols.append('none')
            hatches.append('///')
            edges.append(COL['open'])
        else:
            heights.append(dg)
            if grp == 'ppk2':
                cols.append(COL['ppk2'])
            elif dg > 15:
                cols.append(COL['gs_strong'])
            else:
                cols.append(COL['gs_weak'])
            hatches.append('')
            edges.append('none')

    bars = ax.bar(x, heights, width=0.62, color=cols, edgecolor=edges,
                  hatch=None, linewidth=1.4, zorder=3)
    for b, h in zip(bars, hatches):
        if h:
            b.set_hatch(h)

    for xi, (lab, dg, area, grp, conv) in zip(x, BARS):
        if conv:
            ax.text(xi, dg + 0.45, f'{dg:.1f}', ha='center', fontsize=10,
                    zorder=4)
        else:
            ax.text(xi, 14.4, 'did not\nconverge', ha='center', fontsize=8.5,
                    color='#5F5E5A', zorder=4, style='italic')

    # buried area as a second label line
    ax.set_xticks(x)
    ax.set_xticklabels([f'{lab}\n{a:,} Å²' for (lab, _, a, _, _), _
                        in zip(BARS, x)], fontsize=9.5)

    ax.axvline(3.5, color='#2A2925', lw=0.9, ls=':', zorder=1)
    ax.text(1.5, 23.2, 'GS dodecamer', ha='center', fontsize=11)
    ax.text(5.0, 23.2, 'PPK2 tetramer', ha='center', fontsize=11)

    ax.set_ylabel('interface ΔG density\n(−REU per 1000 Å²)')
    ax.set_ylim(0, 24.5)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    ax.grid(axis='y', alpha=0.25, lw=0.5, zorder=0)

    ax.legend(handles=[
        Patch(facecolor=COL['gs_strong'], label='strongest interface'),
        Patch(facecolor=COL['gs_weak'], label='opens first'),
        Patch(facecolor=COL['ppk2'], label='PPK2, no weak link'),
        Patch(facecolor='none', edgecolor=COL['open'], hatch='///',
              label='not converged'),
    ], frameon=False, fontsize=8.5, ncol=4,
       loc='upper center', bbox_to_anchor=(0.5, -0.17))


def panel_b(ax):
    keys = list(STITCH)
    vals = [STITCH[k] for k in keys]
    x = np.arange(len(keys))

    ax.bar(x, vals, width=0.5, color=[COL['gs_weak'], COL['gs_strong']],
           zorder=3)
    for xi, v in zip(x, vals):
        ax.text(xi, v + 0.7, str(v), ha='center', fontsize=11)

    ax.set_xticks(x)
    ax.set_xticklabels(['within a ring\n(weaker joints)',
                        'across the plane\n(strongest joint)'], fontsize=9.5)
    ax.set_ylabel('genipin-reachable stitches\n(5–16 Å, cross-subunit)')
    ax.set_ylim(0, 33)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    ax.grid(axis='y', alpha=0.25, lw=0.5, zorder=0)

    ax.annotate('', xy=(1, 28.2), xytext=(0, 28.2),
                arrowprops=dict(arrowstyle='<->', color='#5F5E5A', lw=1.1))
    ax.text(0.5, 29.0, 'stitch count runs opposite to interface strength',
            ha='center', fontsize=8.2, color='#5F5E5A')


def main(split=False):
    if split:
        for name, fn, size in (('interface_energetics', panel_a, (7.6, 4.8)),
                               ('interface_stitch_reversal', panel_b, (4.6, 4.4))):
            fig, ax = plt.subplots(figsize=size)
            fn(ax)
            fig.tight_layout()
            fig.savefig(name + '.png', dpi=300, bbox_inches='tight')
            print('wrote ' + name + '.png')
        return

    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=(12.4, 5.0), gridspec_kw={'width_ratios': [2.1, 1]})
    panel_a(ax1)
    panel_b(ax2)
    ax1.set_title('(a)  interface strength', loc='left', fontsize=11)
    ax2.set_title('(b)  crosslinker reach', loc='left', fontsize=11)

    fig.text(0.5, -0.12,
             'ΔG is negative; the axis is negated so taller means stronger. '
             'All converged runs use the fixed-backbone protocol, 5 FastRelax '
             'cycles.\nPPK2 C–D converges backbone-free at −26.7 REU per 1000 Å², '
             'which is a different protocol and is not plotted. Buried area is '
             'the total across both sides.',
             ha='center', fontsize=8.2, color='#5F5E5A')
    fig.tight_layout()
    fig.savefig('interface_energetics.png', dpi=300, bbox_inches='tight')
    print('wrote interface_energetics.png')


if __name__ == '__main__':
    main(split='--split' in sys.argv)
