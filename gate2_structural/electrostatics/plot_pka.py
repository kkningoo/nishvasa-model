import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def read_curve(f):
    ph, dg = [], []
    on = False
    for line in open(f):
        s = line.strip()
        if s.startswith('Free energy of'):
            on = True; continue
        if s.startswith('The pH of optimum'):
            on = False; continue
        if on:
            p = s.split()
            if len(p) == 2:
                try:
                    ph.append(float(p[0])); dg.append(float(p[1]))
                except ValueError:
                    pass
    return ph, dg

def n_titratable(f):
    n, on = 0, False
    for line in open(f):
        s = line.strip()
        if s.startswith('SUMMARY OF THIS PREDICTION'):
            on = True; continue
        if s.startswith('Free energy of'):
            break
        if on:
            p = s.split()
            if len(p) >= 4 and p[0] in ('ASP','GLU','HIS','CYS','TYR','LYS','ARG','N+','C-'):
                n += 1
    return n

files = [('ppk2_prot.pka', 'PPK2 tetramer',  '#1D9E75'),
         ('gs_tmp.pka',    'GS dodecamer',   '#534AB7')]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2))

for f, label, col in files:
    ph, dg = read_curve(f)
    n = n_titratable(f)
    lo = min(dg)
    rel = [(d - lo) / n for d in dg]
    ax1.plot(ph, dg, color=col, lw=1.8, label=f'{label} (n={n})')
    ax2.plot(ph, rel, color=col, lw=1.8, label=label)
    imin = dg.index(lo)
    ax2.plot(ph[imin], 0, 'o', color=col, ms=5)

for ax in (ax1, ax2):
    ax.set_xlabel('pH')
    ax.set_xlim(2, 12)
    ax.axvspan(6.5, 8.0, color='#888780', alpha=0.12, lw=0)
    ax.grid(alpha=0.25, lw=0.5)
    ax.legend(frameon=False, fontsize=9)

ax1.set_ylabel('PROPKA $\\Delta G$ (kcal/mol)')
ax1.set_title('Raw output', fontsize=10)
ax2.set_ylabel('kcal/mol per titratable residue')
ax2.set_ylim(-0.01, 0.30)
ax2.set_title('Relative to each protein\'s own optimum', fontsize=10)

fig.text(0.5, -0.02,
         'Shaded band = cartridge operating range. Metals included in the GS calculation.',
         ha='center', fontsize=8, color='#5F5E5A')
fig.tight_layout()
fig.savefig('pka_stability.png', dpi=300, bbox_inches='tight')
print('wrote pka_stability.png')
