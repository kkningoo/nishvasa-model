"""
Model 5 — Honeycomb Scaffold Cell Geometry
Picks the hexagonal cell inradius r by maximising VRE (volumetric
removal efficiency = GSA × ε_open). An illustrative Design Score
heuristic is shown but NOT used to set r.

All lengths in m except where _mm suffix is used.
"""

import numpy as np
import matplotlib.pyplot as plt

# --- Parameters ---

L_COAT = 24e-6     # m
T_WALL = 0.8e-3    # m

# analytical VRE optimum
R_STAR = 4 * L_COAT + 1.5 * T_WALL
print(f"r* = 4·L_coat + 1.5·t_w = {R_STAR*1e3:.3f} mm")

# algebraic check: at r*, r_open/r_out = 3/4 exactly, so eps_open = 9/16 = 56.25%
r_open_star = R_STAR - L_COAT
r_out_star = R_STAR + T_WALL / 2
ratio_check = r_open_star / r_out_star
eps_check = ratio_check**2
print(f"At r*: r_open/r_out = {ratio_check:.4f} (expect 0.75), "
      f"ε_open = {eps_check:.4f} (expect 0.5625 = 9/16)")


def cell_geometry(r):
    """Returns r_open, r_out, A_open, A_cell, P, D_h, eps_open, GSA, VRE."""
    r_open = r - L_COAT
    r_out = r + T_WALL / 2
    A_open = 2 * np.sqrt(3) * r_open**2
    A_cell = 2 * np.sqrt(3) * r_out**2
    P = 4 * np.sqrt(3) * r_open
    D_h = 2 * r_open
    eps_open = A_open / A_cell
    GSA = P / A_cell
    VRE = GSA * eps_open  # = 2*r_open^3 / r_out^4
    return dict(r_open=r_open, r_out=r_out, A_open=A_open, A_cell=A_cell,
                P=P, D_h=D_h, eps_open=eps_open, GSA=GSA, VRE=VRE)


# --- Sweep r from 0.3 mm to 5 mm ---

r_sweep = np.linspace(0.3e-3, 5e-3, 2000)
r_mm = r_sweep * 1e3
vre_arr = np.array([cell_geometry(r)["VRE"] for r in r_sweep])

# normalise VRE to peak
vre_peak = vre_arr.max()
vre_norm = vre_arr / vre_peak
i_vrepeak = np.argmax(vre_arr)
r_vrepeak_mm = r_mm[i_vrepeak]

# 97% plateau bounds
above_97 = r_mm[vre_norm >= 0.97]
plateau_lo = above_97[0]
plateau_hi = above_97[-1]

# --- Design Score (illustrative heuristic) ---

R_ACCESS = 0.45e-3   # m — judgement value
D_MIN = 0.8e-3       # m — FDM minimum feature

F_coat = 1.0 - np.exp(-r_sweep / R_ACCESS)
F_fdm = 1.0 - (D_MIN / (2 * r_sweep))**2
design_score = vre_arr * F_coat * F_fdm
ds_peak = design_score.max()
ds_norm = design_score / ds_peak
i_dspeak = np.argmax(design_score)
r_dspeak_mm = r_mm[i_dspeak]

# --- Print table ---

print(f"\nVRE peak at r = {r_vrepeak_mm:.3f} mm (= r*)")
print(f"Design Score peak at r = {r_dspeak_mm:.2f} mm")
print(f"97% VRE plateau: {plateau_lo:.1f} to {plateau_hi:.1f} mm\n")

test_radii_mm = [0.5, 1.0, 1.5, 2.0, 3.0]
print(f"{'r (mm)':>8}  {'ε_open (%)':>10}  {'D_h (mm)':>8}  {'VRE (rel)':>10}  {'DS (rel)':>10}")
for rm in test_radii_mm:
    g = cell_geometry(rm * 1e-3)
    vn = g["VRE"] / vre_peak
    dsn = (g["VRE"] * (1.0 - np.exp(-rm*1e-3 / R_ACCESS)) *
           (1.0 - (D_MIN / (2 * rm * 1e-3))**2)) / ds_peak
    print(f"{rm:8.1f}  {g['eps_open']*100:10.1f}  {g['D_h']*1e3:8.2f}  {vn:10.2f}  {dsn:10.2f}")

# selected r = 1.5 mm
g_sel = cell_geometry(1.5e-3)
vn_sel = g_sel["VRE"] / vre_peak
print(f"\nSelected r = 1.5 mm retains {vn_sel*100:.1f}% of VRE peak.")
print(f"Geometry: r_open = {g_sel['r_open']*1e3:.3f} mm, D_h = {g_sel['D_h']*1e3:.2f} mm, "
      f"A_open = {g_sel['A_open']*1e6:.2f} mm², A_cell = {g_sel['A_cell']*1e6:.1f} mm², "
      f"P = {g_sel['P']*1e3:.1f} mm, ε_open = {g_sel['eps_open']*100:.1f}%")

# circumscribed-vs-hex area error check
hex_to_circle = 2 * np.pi / (3 * np.sqrt(3))
print(f"\nCircumscribed-circle / hex area ratio = {hex_to_circle:.3f} "
      f"(a {(hex_to_circle - 1)*100:.0f}% overestimate)")


# ============================================================
# Figure: normalised VRE and Design Score vs r
# ============================================================

fig, ax = plt.subplots(figsize=(9, 5.5), constrained_layout=True)

# 97% plateau shading
ax.axhspan(0.97, 1.0, color="#0072B2", alpha=0.08, zorder=0)
ax.axhline(0.97, color="#0072B2", linestyle=":", alpha=0.3)
ax.text(4.5, 0.975, "≥ 97% of VRE peak", color="#0072B2", fontsize=8, alpha=0.6,
        ha="right", va="bottom")

ax.plot(r_mm, vre_norm, color="#0072B2", linewidth=2.5, label="VRE (normalised)")
ax.plot(r_mm, ds_norm, color="#E69F00", linewidth=2, linestyle="--",
        label="Design Score (illustrative)")

# r* mark
ax.axvline(r_vrepeak_mm, color="#0072B2", linestyle=":", alpha=0.5)
ax.annotate(f"VRE peak\nr* = {r_vrepeak_mm:.2f} mm",
            xy=(r_vrepeak_mm, 1.0), xytext=(r_vrepeak_mm + 0.5, 0.88),
            arrowprops=dict(arrowstyle="->", color="#0072B2"),
            fontsize=9, color="#0072B2")

# selected r = 1.5 mm
ax.plot(1.5, vn_sel, 'o', color="#0072B2", markersize=9, zorder=5)
ax.annotate(f"Selected\nr = 1.5 mm\n({vn_sel*100:.1f}% of peak)",
            xy=(1.5, vn_sel), xytext=(2.2, 0.80),
            arrowprops=dict(arrowstyle="->", color="#0072B2"),
            fontsize=9, color="#0072B2")

# DS peak
ax.annotate(f"DS peak\n{r_dspeak_mm:.2f} mm",
            xy=(r_dspeak_mm, 1.0), xytext=(r_dspeak_mm + 0.6, 0.72),
            arrowprops=dict(arrowstyle="->", color="#E69F00"),
            fontsize=9, color="#E69F00")

ax.set_xlabel("Cell inradius r (mm)", fontsize=11)
ax.set_ylabel("Normalised value", fontsize=11)
ax.set_title("VRE and Design Score vs cell inradius\n"
             f"(L_coat = {L_COAT*1e6:.0f} µm, t_w = {T_WALL*1e3:.1f} mm)", fontsize=12)
ax.legend(frameon=True, fontsize=10, loc="lower right")
ax.set_xlim(0.3, 5.0)
ax.set_ylim(0, 1.08)

fig.savefig("model5_scaffold_geometry.png", dpi=200, bbox_inches="tight")
plt.close(fig)

print("\nConclusion: r = 1.5 mm is inside the ≥ 97% VRE plateau (1.0–1.7 mm). "
      "It is selected for its higher open fraction and printability, not at either peak.")
