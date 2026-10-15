"""Two figures for GS variant selection, to spec.

Fig 1  effective rate against accumulated glutamine, log x to 1000 mM,
       all six variants, run range shaded
Fig 2  heatmap of v_H195Y / v_M68I over ATP x glutamine, showing that no
       crossover exists anywhere in the swept domain

Everything derives from Fisher & Wray (2006) Tables 3 and 4 plus the
operating conditions. Nothing is entered by hand downstream of those.

    python3 variant_figures2.py
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm

# ---- inputs: Fisher & Wray (2006), Tables 3 and 4 ------------------------
#          Km_Glu  Km_ATP  a      IC50_Gln  bounded?
VARIANTS = {
    'Wild type': (27.0,  2.4,  2.7,   2.4,  False),
    'E65G':      (7.3,   3.6,  0.35,  140,  True),
    'S66P':      (130.0, 23.0, 2.4,   140,  True),
    'M68I':      (56.0,  4.1,  2.4,   139,  False),
    'H195Y':     (30.0,  7.0,  4.7,   115,  False),
    'P318S':     (3.3,   0.81, 0.86,  140,  True),
}

GLU_ASSAY = 150.0
GLU_OP    = 150.0
ATP_OP    = 3.8
NH3_OP    = 2.95
KM_NH3    = 0.4

# Glutamine traversed during operation. Set from the run you want to show:
# a single batch accumulates far less than a multi-day reservoir cycle.
GLN_RUN = (0.0, 20.0)      # [FILL: confirm against the reservoir model]

# Sweep domain for Fig 2
ATP_RANGE = (0.1, 50.0)
GLN_RANGE = (0.0, 1000.0)

COL = {
    'H195Y':     '#534AB7',
    'M68I':      '#1D9E75',
    'Wild type': '#C0392B',
    'P318S':     '#D85A30',
    'S66P':      '#A8A59E',
    'E65G':      '#C9C6BE',
}
ORDER = ['H195Y', 'M68I', 'P318S', 'Wild type', 'S66P', 'E65G']


def ki(km_glu, ic50):
    return ic50 / (1.0 + GLU_ASSAY / km_glu)


def rate(km_glu, km_atp, a, k_i, gln, atp=ATP_OP):
    f_atp = atp / (km_atp + atp)
    f_nh3 = NH3_OP / (KM_NH3 + NH3_OP)
    f_glu = GLU_OP / (km_glu * (1.0 + gln / k_i) + GLU_OP)
    return a * f_atp * f_nh3 * f_glu


def figure1():
    """Log x so the wild-type collapse is visible rather than a vertical drop."""
    gln = np.logspace(-1, 3, 800)          # 0.1 to 1000 mM
    fig, ax = plt.subplots(figsize=(7.4, 4.8))

    lo, hi = GLN_RUN
    ax.axvspan(max(lo, 0.1), hi, color='#534AB7', alpha=0.07, lw=0, zorder=1)
    ax.text(np.sqrt(max(lo, 0.1) * hi), 1.30, 'range traversed\nduring a run',
            fontsize=8, color='#5F5E5A', ha='center', va='top')

    for name in ORDER:
        kg, ka, a, ic50, _ = VARIANTS[name]
        k = ki(kg, ic50)
        y = rate(kg, ka, a, k, gln)
        lead = name in ('H195Y', 'M68I', 'Wild type')
        ax.plot(gln, y, color=COL[name], lw=2.2 if lead else 1.3,
                label=name, zorder=4 if name == 'H195Y' else 3,
                alpha=1.0 if lead else 0.85)

    ax.set_xscale('log')
    ax.set_xlim(0.1, 1000)
    ax.set_ylim(0, 1.38)
    ax.set_xlabel('accumulated glutamine (mM, log scale)')
    ax.set_ylabel('effective specific activity\n'
                  '($\\mu$mol min$^{-1}$ mg$^{-1}$)')
    ax.legend(frameon=False, fontsize=9, ncol=1, loc='upper right')
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    ax.grid(alpha=0.25, lw=0.5)

    fig.text(0.5, -0.03,
             f'[Glu] = {GLU_OP:.0f} mM, [ATP] = {ATP_OP} mM, '
             f'[NH$_3$] = {NH3_OP} mM, 25 $\\degree$C. Literature kinetics on '
             f'free enzyme; no error bars available.',
             ha='center', fontsize=8, color='#5F5E5A')
    fig.tight_layout()
    fig.savefig('fig1_rate_vs_glutamine.png', dpi=300, bbox_inches='tight')
    print('wrote fig1_rate_vs_glutamine.png')


def figure2():
    """Ratio over the full operating domain. No contour at 1 means no crossover."""
    kH = ki(*[VARIANTS['H195Y'][i] for i in (0, 3)])
    kM = ki(*[VARIANTS['M68I'][i] for i in (0, 3)])
    hH = VARIANTS['H195Y']
    hM = VARIANTS['M68I']

    atp = np.logspace(np.log10(ATP_RANGE[0]), np.log10(ATP_RANGE[1]), 500)
    gln = np.linspace(*GLN_RANGE, 500)
    A, G = np.meshgrid(atp, gln)
    R = (rate(hH[0], hH[1], hH[2], kH, G, A) /
         rate(hM[0], hM[1], hM[2], kM, G, A))

    i, j = np.unravel_index(R.argmin(), R.shape)
    rmin, amin, gmin = R[i, j], A[i, j], G[i, j]
    rop = rate(hH[0], hH[1], hH[2], kH, GLN_RUN[1], ATP_OP) / \
          rate(hM[0], hM[1], hM[2], kM, GLN_RUN[1], ATP_OP)

    fig, ax = plt.subplots(figsize=(7.0, 4.8))
    im = ax.pcolormesh(A, G, R, cmap='PuBuGn', shading='auto',
                       norm=LogNorm(vmin=1.0, vmax=R.max()))

    cs = ax.contour(A, G, R, levels=[1.2, 1.4, 1.6, 1.8, 2.0],
                    colors='white', linewidths=0.8, alpha=0.8)
    ax.clabel(cs, inline=True, fontsize=7.5, fmt='%.1f')

    # the ratio = 1 contour, which does not exist in this domain
    ax.contour(A, G, R, levels=[1.0], colors='#C0392B', linewidths=2.0)

    ax.scatter([ATP_OP], [GLN_RUN[1]], s=70, color='#C0392B', zorder=5,
               edgecolor='white', lw=1.3)
    ax.annotate(f'operating point\n{rop:.2f}', (ATP_OP, GLN_RUN[1]),
                xytext=(14, 22), textcoords='offset points',
                fontsize=8.5, color='#C0392B',
                arrowprops=dict(arrowstyle='-', color='#C0392B', lw=0.9))

    ax.scatter([amin], [gmin], s=70, marker='s', color='#2A2925', zorder=5,
               edgecolor='white', lw=1.3)
    ax.annotate(f'minimum\n{rmin:.2f}', (amin, gmin),
                xytext=(16, -28), textcoords='offset points',
                fontsize=8.5, color='#2A2925',
                arrowprops=dict(arrowstyle='-', color='#2A2925', lw=0.9))

    ax.set_xscale('log')
    ax.set_xlabel('[ATP] (mM, log scale)')
    ax.set_ylabel('accumulated glutamine (mM)')
    cb = fig.colorbar(im, ax=ax, pad=0.02)
    cb.set_label('$v_{\\mathrm{H195Y}}\\,/\\,v_{\\mathrm{M68I}}$')
    cb.ax.axhline(1.0, color='#C0392B', lw=2.0)

    fig.text(0.5, -0.03,
             f'Ratio ranges {rmin:.2f} to {R.max():.2f} across the domain. The '
             f'crossover contour at 1.0 (red) does not appear:\nno combination '
             f'of substrate and product concentration exists at which M68I '
             f'overtakes H195Y.',
             ha='center', fontsize=8, color='#5F5E5A')
    fig.tight_layout()
    fig.savefig('fig2_ratio_heatmap.png', dpi=300, bbox_inches='tight')
    print('wrote fig2_ratio_heatmap.png')
    print(f'  minimum {rmin:.4f} at ATP {amin:.2f} mM, Gln {gmin:.0f} mM')
    print(f'  maximum {R.max():.3f}')
    print(f'  operating point ({ATP_OP} mM, {GLN_RUN[1]:.0f} mM): {rop:.3f}')
    print(f'  fraction of domain below 1.0: {(R < 1).mean():.4f}')


if __name__ == '__main__':
    figure1()
    figure2()
