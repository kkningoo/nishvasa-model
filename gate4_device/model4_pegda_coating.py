"""
Model 4 — PEGDA Coating Thickness
Trade-off between enzyme utilisation (η_coat) and reactive volume
fraction (ε_coat) in the hydrogel coating on scaffold walls.

Sweep in the real 1.5 mm scaffold cell. J = η_coat × ε_coat is shown
but NOT used to justify the 24 µm selection — that is retained on
enzyme-utilisation grounds (η_coat = 0.985).

All lengths in m except where _um or _mm suffix is used.
"""

import numpy as np
import matplotlib.pyplot as plt

# --- Parameters (same ATP basis as Model 3) ---

A_GS = 0.0783333             # mol s⁻¹ kg⁻¹
M_LOAD = 2.0                 # kg/m³
KM_PLUS_C_ATP = 10.5         # mM = mol/m³
D_ATP_GEL = 1.85e-10         # m²/s

R_CELL = 1.5e-3              # m, hex cell inradius

# φ_coat = L_coat × k, where k = √(a·M / ((Km+c)·D))
k = np.sqrt(A_GS * M_LOAD / (KM_PLUS_C_ATP * D_ATP_GEL))
saturation_depth = 1.0 / k   # asymptote of η_coat × L_coat

print(f"Rate constant k = {k:.2f} /m = {k*1e-3:.4f} /mm")
print(f"Saturation depth 1/k = {saturation_depth*1e6:.0f} µm")


def eta_coat(L_coat):
    phi = L_coat * k
    return np.where(phi > 1e-6, np.tanh(phi) / phi, 1.0)


def eps_coat(L_coat, r):
    return 1.0 - (1.0 - L_coat / r)**2


def J_metric(L_coat, r):
    return eta_coat(L_coat) * eps_coat(L_coat, r)


# --- Sweep ---

L_sweep = np.linspace(1e-6, 700e-6, 5000)  # 1–700 µm

eta_arr = eta_coat(L_sweep)
eps_arr = eps_coat(L_sweep, R_CELL)
J_arr = J_metric(L_sweep, R_CELL)
active_depth_arr = eta_arr * L_sweep

# selected point
L_sel = 24e-6
eta_sel = eta_coat(np.array([L_sel]))[0]
eps_sel = eps_coat(np.array([L_sel]), R_CELL)[0]
J_sel = J_metric(np.array([L_sel]), R_CELL)[0]

# J peak
i_peak = np.argmax(J_arr)
L_peak = L_sweep[i_peak]
J_peak_val = J_arr[i_peak]
eta_peak = eta_arr[i_peak]
eps_peak = eps_arr[i_peak]

print(f"\nAt L_coat = 24 µm:")
print(f"  η_coat  = {eta_sel:.4f}")
print(f"  ε_coat  = {eps_sel:.4f} ({eps_sel*100:.1f}%)")
print(f"  J       = {J_sel:.4f}")
print(f"  J / J_peak = {J_sel/J_peak_val*100:.1f}%")

print(f"\nJ peak at L_coat = {L_peak*1e6:.0f} µm:")
print(f"  η_coat = {eta_peak:.2f}, ε_coat = {eps_peak:.2f}, J = {J_peak_val:.3f}")

# check with open radius
R_OPEN = R_CELL - L_sel
eps_sel_open = eps_coat(np.array([L_sel]), R_OPEN)[0]
J_sel_open = eta_sel * eps_sel_open
print(f"\nUsing r_open = {R_OPEN*1e3:.3f} mm instead of r = {R_CELL*1e3:.1f} mm:")
print(f"  ε_coat = {eps_sel_open:.4f}, J = {J_sel_open:.4f}")
print(f"  J change = {abs(J_sel_open - J_sel)/J_sel*100:.1f}%")

print(f"\nThe coating thickness is retained at 24 µm on enzyme-utilisation grounds "
      f"(η_coat = {eta_sel:.3f}, i.e. {(1-eta_sel)*100:.1f}% activity cost), "
      f"NOT because J is maximised there.")


# ============================================================
# Figure: two panels
# ============================================================

fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(14, 5.5), constrained_layout=True)

L_um = L_sweep * 1e6

# --- Panel (a): J and ε_coat ---

color_j = "#0072B2"
color_eps = "#D55E00"

ax_a.plot(L_um, J_arr, color=color_j, linewidth=2, label=r"$J = \eta_{\mathrm{coat}} \times \varepsilon_{\mathrm{coat}}$")
ax_a_twin = ax_a.twinx()
ax_a_twin.plot(L_um, eps_arr * 100, color=color_eps, linewidth=2, linestyle="--",
               label=r"$\varepsilon_{\mathrm{coat}}$")

# peak mark
ax_a.plot(L_peak * 1e6, J_peak_val, 'o', color=color_j, markersize=8, zorder=5)
ax_a.annotate(f"J peak at {L_peak*1e6:.0f} µm\n(J = {J_peak_val:.3f})",
              xy=(L_peak * 1e6, J_peak_val),
              xytext=(L_peak * 1e6 + 120, J_peak_val - 0.02),
              arrowprops=dict(arrowstyle="->", color=color_j),
              fontsize=9, color=color_j)

# 24 µm mark
ax_a.plot(24, J_sel, 's', color=color_j, markersize=8, zorder=5)
ax_a.annotate(f"24 µm\n(J = {J_sel:.3f})",
              xy=(24, J_sel),
              xytext=(80, J_sel + 0.02),
              arrowprops=dict(arrowstyle="->", color=color_j),
              fontsize=9, color=color_j)

ax_a.set_xlabel(r"Coating thickness $L_{\mathrm{coat}}$ (µm)", fontsize=11)
ax_a.set_ylabel(r"$J$", fontsize=11, color=color_j)
ax_a_twin.set_ylabel(r"$\varepsilon_{\mathrm{coat}}$ (%)", fontsize=11, color=color_eps)
ax_a.set_title(r"(a) $J$ and $\varepsilon_{\mathrm{coat}}$ vs coating thickness", fontsize=12)
ax_a.tick_params(axis="y", colors=color_j)
ax_a_twin.tick_params(axis="y", colors=color_eps)
ax_a.set_xlim(0, 700)
ax_a.set_ylim(0, 0.16)
ax_a_twin.set_ylim(0, 100)

# combined legend
lines_a = ax_a.get_lines() + ax_a_twin.get_lines()
labels_a = [l.get_label() for l in lines_a]
ax_a.legend(lines_a, labels_a, frameon=True, fontsize=9, loc="center right")

# --- Panel (b): η_coat and active depth ---

color_eta = "#0072B2"
color_depth = "#009E73"

ax_b.plot(L_um, eta_arr, color=color_eta, linewidth=2, label=r"$\eta_{\mathrm{coat}}$")
ax_b_twin = ax_b.twinx()
ax_b_twin.plot(L_um, active_depth_arr * 1e6, color=color_depth, linewidth=2, linestyle="--",
               label=r"Active depth $\eta_{\mathrm{coat}} L$")

# 1/k asymptote
ax_b_twin.axhline(saturation_depth * 1e6, color=color_depth, linestyle=":", alpha=0.6)
ax_b_twin.text(600, saturation_depth * 1e6 + 3, f"1/k = {saturation_depth*1e6:.0f} µm",
               color=color_depth, fontsize=9, ha="right", va="bottom")

ax_b.set_xlabel(r"Coating thickness $L_{\mathrm{coat}}$ (µm)", fontsize=11)
ax_b.set_ylabel(r"$\eta_{\mathrm{coat}}$", fontsize=11, color=color_eta)
ax_b_twin.set_ylabel(r"Active depth $\eta L$ (µm)", fontsize=11, color=color_depth)
ax_b.set_title(r"(b) $\eta_{\mathrm{coat}}$ and active depth", fontsize=12)
ax_b.tick_params(axis="y", colors=color_eta)
ax_b_twin.tick_params(axis="y", colors=color_depth)
ax_b.set_xlim(0, 700)
ax_b.set_ylim(0, 1.05)
ax_b_twin.set_ylim(0, 130)

lines_b = ax_b.get_lines() + ax_b_twin.get_lines()
labels_b = [l.get_label() for l in lines_b]
ax_b.legend(lines_b, labels_b, frameon=True, fontsize=9, loc="center right")

fig.savefig("model4_pegda_coating.png", dpi=200, bbox_inches="tight")
plt.close(fig)
