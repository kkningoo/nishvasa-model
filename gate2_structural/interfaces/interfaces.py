import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

# GS, fixed-backbone protocol
gs = [("A-H\ncross-ring", 21.20, "#534AB7"),
      ("A-F\nintra-ring", 12.85, "#888780"),
      ("A-B\nintra-ring", 10.28, "#888780")]

# PPK2, fixed-backbone protocol (C-D omitted - required backbone-free relax)
p2 = [("A-C", 19.96, "#1D9E75"),
      ("A-B", 18.34, "#1D9E75")]

fig, ax = plt.subplots(figsize=(7, 4.2))
labels = [x[0] for x in gs+p2]
vals   = [x[1] for x in gs+p2]
cols   = [x[2] for x in gs+p2]
xs = np.arange(len(vals))
ax.bar(xs, vals, color=cols, width=0.62)
ax.set_xticks(xs); ax.set_xticklabels(labels, fontsize=9)
ax.set_ylabel("interface $\\Delta$G density\n(-REU per 1000 $\\AA^2$)")
ax.set_ylim(0, 24)
for x, v in zip(xs, vals):
    ax.text(x, v+0.4, f"{v:.1f}", ha="center", fontsize=9)
ax.axvline(2.5, color="black", lw=0.8, ls=":")
ax.text(1.0, 23, "GS dodecamer", ha="center", fontsize=10)
ax.text(3.5, 23, "PPK2 tetramer", ha="center", fontsize=10)
for s in ("top","right"): ax.spines[s].set_visible(False)
ax.grid(axis="y", alpha=0.25, lw=0.5)
fig.tight_layout()
fig.savefig("interface_energetics.png", dpi=300)
print("wrote interface_energetics.png")
