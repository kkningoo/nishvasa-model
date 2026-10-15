"""
Model 3 — Bead Effectiveness Factor
Checks whether intraparticle diffusion limits the immobilized enzyme at
the corrected bead radius (10 µm, constrained by coating geometry, not
by diffusion). Compares sphere/cylinder/slab shapes.

All lengths in m, diffusivities in m²/s, concentrations in mM where
needed for Km+c, rates in mol/(s·m³).
"""

import numpy as np
from scipy.special import iv as bessel_iv  # modified Bessel I_v
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

# --- Physical parameters ---

A_GS = 0.0783333             # mol s⁻¹ kg⁻¹  (4.7 µmol/min/mg)
M_LOAD = 2.0                 # kg/m³
A_M = A_GS * M_LOAD          # mol s⁻¹ m⁻³ = intrinsic volumetric rate constant

D_FREE_ATP = 3.7e-10         # m²/s
D_FREE_NH3 = 1.6e-9
D_FREE_GLU = 7.6e-10         # approximate, for glutamate in water
GEL_RETARDATION = 0.5

D_ATP_GEL = GEL_RETARDATION * D_FREE_ATP   # 1.85e-10
D_NH3_GEL = GEL_RETARDATION * D_FREE_NH3   # 8.0e-10
D_GLU_GEL = GEL_RETARDATION * D_FREE_GLU   # 3.8e-10

# Km + c for each substrate at operating conditions (mM -> mol/m³ for Thiele, but
# ratio a*M / ((Km+c)*D) has Km+c in same units as a*M*L² in denominator —
# actually Km+c in mol/m³ and a*M in mol/(s·m³), D in m²/s, L in m)
# Convert mM to mol/m³: × 1.0

KM_PLUS_C = {
    "ATP":       (7.0 + 3.5, D_ATP_GEL),      # mM = mol/m³
    "NH3":       (0.4 + 2.95, 8.2e-10),
    "Glutamate": (30.0 + 150.0, 3.8e-10),
}

L_COAT = 24e-6               # m
L_BEAD_SELECTED = 10e-6      # m


def thiele_modulus(L, a_M, km_plus_c, D):
    """φ = L √(a·M / ((Km+c)·D)), with Km+c in mol/m³."""
    return L * np.sqrt(a_M / (km_plus_c * D))


def eta_sphere(phi):
    """Closed-form effectiveness factor for a sphere, with series for small φ."""
    phi = np.atleast_1d(np.float64(phi))
    eta = np.ones_like(phi)
    large = phi > 1e-2
    eta[large] = (3.0 / phi[large]) * (1.0 / np.tanh(phi[large]) - 1.0 / phi[large])
    small = ~large
    eta[small] = 1.0 - phi[small]**2 / 15.0
    return eta


def eta_cylinder(phi):
    phi = np.atleast_1d(np.float64(phi))
    eta = np.ones_like(phi)
    large = phi > 1e-6
    eta[large] = (2.0 / phi[large]) * bessel_iv(1, phi[large]) / bessel_iv(0, phi[large])
    return eta


def eta_slab(phi):
    phi = np.atleast_1d(np.float64(phi))
    eta = np.ones_like(phi)
    large = phi > 1e-6
    eta[large] = np.tanh(phi[large]) / phi[large]
    return eta


# --- Print Thiele modulus for each substrate at L = 10 µm ---

print("Thiele modulus at L = 10 µm, M = 2.0 kg/m³:")
for substrate, (kmc, D) in KM_PLUS_C.items():
    phi = thiele_modulus(L_BEAD_SELECTED, A_M, kmc, D)
    eta = eta_sphere(np.array([phi]))[0]
    print(f"  {substrate:12s}: Km+c = {kmc:6.2f} mM, D = {D:.2e} m²/s, "
          f"φ = {phi:.4f}, η_sphere = {eta:.4f}")

# --- Diffusivity audit ---

# The old tuned constant implied these diffusivities (back-calculated):
D_implied_ATP = 4.8e-12
D_implied_NH3 = 1.5e-11
print(f"\nDiffusivity audit (withdrawn tuned constant):")
print(f"  Implied D_ATP = {D_implied_ATP:.1e} m²/s  vs  free = {D_FREE_ATP:.1e} "
      f"({D_FREE_ATP/D_implied_ATP:.0f}× higher)")
print(f"  Implied D_NH3 = {D_implied_NH3:.1e} m²/s  vs  free = {D_FREE_NH3:.1e} "
      f"({D_FREE_NH3/D_implied_NH3:.0f}× higher)")
print(f"  Both are 40–100× below the plausible gel range. Tuned constant discarded.")

# --- Bead radius constraint ---

L_max = L_COAT / 2
print(f"\nBead radius constraint: L_bead ≤ L_coat/2 = {L_max*1e6:.0f} µm")
print(f"Selected: L_bead = {L_BEAD_SELECTED*1e6:.0f} µm (leaves margin for coating variation)")
print(f"Leaching caveat: 10 µm is below the 34–500 µm chitosan microsphere range. "
      f"Enzyme retention at this size is a wet-lab question, not captured by η.")

# --- Shape comparison at phi = 5 ---

phi_check = 5.0
print(f"\nShape comparison at φ = {phi_check}:")
print(f"  η_sphere   = {eta_sphere(np.array([phi_check]))[0]:.2f}")
print(f"  η_cylinder = {eta_cylinder(np.array([phi_check]))[0]:.2f}")
print(f"  η_slab     = {eta_slab(np.array([phi_check]))[0]:.2f}")

# η at M = 30 kg/m³, L = 10 µm
phi_30 = thiele_modulus(L_BEAD_SELECTED, A_GS * 30.0, 10.5, D_ATP_GEL)
eta_30 = eta_sphere(np.array([phi_30]))[0]
print(f"\nAt M = 30 kg/m³, L = 10 µm: φ = {phi_30:.3f}, η_sphere = {eta_30:.4f}")


# ============================================================
# Figure 1: η_bead vs assumed ATP diffusivity
# ============================================================

D_sweep = np.logspace(-12, -9, 500)  # m²/s
phi_D = thiele_modulus(L_BEAD_SELECTED, A_M, 10.5, D_sweep)
eta_D = eta_sphere(phi_D)

fig1, ax1 = plt.subplots(figsize=(8, 5), constrained_layout=True)
ax1.semilogx(D_sweep, eta_D, color="#0072B2", linewidth=2)

# plausible gel range: 0.3–0.8 × free
D_gel_lo = 0.3 * D_FREE_ATP
D_gel_hi = 0.8 * D_FREE_ATP
ax1.axvspan(D_gel_lo, D_gel_hi, color="#009E73", alpha=0.15, label="Plausible gel range (0.3–0.8 × free)")

# stated value
ax1.axvline(D_ATP_GEL, color="#009E73", linestyle="--", alpha=0.7, label=f"Stated D = {D_ATP_GEL:.2e}")

# withdrawn value
ax1.axvline(D_implied_ATP, color="#D55E00", linestyle=":", alpha=0.7,
            label=f"Withdrawn (tuned) D = {D_implied_ATP:.1e}")

ax1.set_xlabel(r"Assumed $D_{\mathrm{ATP,gel}}$ (m²/s)", fontsize=11)
ax1.set_ylabel(r"$\eta_{\mathrm{bead}}$", fontsize=11)
ax1.set_title(r"Bead effectiveness factor vs assumed ATP diffusivity"
              "\n" r"($L$ = 10 µm, $M$ = 2.0 kg/m³)", fontsize=12)
ax1.legend(frameon=True, fontsize=9, loc="lower right")
ax1.set_ylim(0.9, 1.005)
fig1.savefig("model3_fig1_eta_vs_diffusivity.png", dpi=200, bbox_inches="tight")
plt.close(fig1)


# ============================================================
# Figure 2: η vs Thiele modulus for three shapes
# ============================================================

phi_range = np.logspace(-2, 1.5, 500)

fig2, ax2 = plt.subplots(figsize=(8, 5), constrained_layout=True)
ax2.semilogx(phi_range, eta_sphere(phi_range), color="#0072B2", linewidth=2, label="Sphere")
ax2.semilogx(phi_range, eta_cylinder(phi_range), color="#E69F00", linewidth=2, label="Cylinder")
ax2.semilogx(phi_range, eta_slab(phi_range), color="#D55E00", linewidth=2, label="Slab")

# operating point
phi_op = thiele_modulus(L_BEAD_SELECTED, A_M, 10.5, D_ATP_GEL)
ax2.axvline(phi_op, color="gray", linestyle="--", alpha=0.5)
ax2.annotate(f"Operating φ = {phi_op:.3f}",
             xy=(phi_op, 0.95), xytext=(phi_op * 5, 0.75),
             arrowprops=dict(arrowstyle="->", color="gray"), fontsize=9, color="gray")

ax2.set_xlabel(r"Thiele modulus $\phi$", fontsize=11)
ax2.set_ylabel(r"$\eta$", fontsize=11)
ax2.set_title("Effectiveness factor vs Thiele modulus (Aris 1975)", fontsize=12)
ax2.legend(frameon=True, fontsize=10)
ax2.set_ylim(0, 1.05)
fig2.savefig("model3_fig2_eta_vs_phi_shapes.png", dpi=200, bbox_inches="tight")
plt.close(fig2)


# ============================================================
# Figure 3: η map over bead radius and loading
# ============================================================

L_range = np.logspace(np.log10(5e-6), np.log10(200e-6), 200)   # 5–200 µm
M_range = np.logspace(np.log10(1.0), np.log10(30.0), 200)       # 1–30 kg/m³

L_grid, M_grid = np.meshgrid(L_range, M_range)
phi_grid = thiele_modulus(L_grid, A_GS * M_grid, 10.5, D_ATP_GEL)
eta_grid = eta_sphere(phi_grid)

fig3, ax3 = plt.subplots(figsize=(8, 6), constrained_layout=True)
pcm = ax3.pcolormesh(L_range * 1e6, M_range, eta_grid, cmap="viridis", shading="auto",
                     vmin=0.9, vmax=1.0)
cbar = fig3.colorbar(pcm, ax=ax3, label=r"$\eta_{\mathrm{sphere}}$")

# coating-fit limit at 12 µm
ax3.axvline(12.0, color="white", linestyle="--", linewidth=1.5)
ax3.text(12.5, 28, "L ≤ 12 µm\n(coating fit)", color="white", fontsize=8, va="top")

# operating point
ax3.plot(10, 2.0, '*', color="white", markersize=14, markeredgecolor="black", markeredgewidth=0.5)
ax3.text(13, 2.0, "Selected\n(10 µm, 2 kg/m³)", color="white", fontsize=8, va="center")

ax3.set_xlabel("Bead radius (µm)", fontsize=11)
ax3.set_ylabel(r"Enzyme loading $M$ (kg/m³)", fontsize=11)
ax3.set_title(r"$\eta_{\mathrm{sphere}}$ over bead radius and loading (ATP-limited)", fontsize=12)
ax3.set_xscale("log")
ax3.set_yscale("log")

fig3.savefig("model3_fig3_eta_map.png", dpi=200, bbox_inches="tight")
plt.close(fig3)

print("\nConclusion: At L = 10 µm, M = 2.0 kg/m³, η_bead = 0.9995. Diffusion does not "
      "limit the bead. Bead size is set by coating geometry (L ≤ 12 µm), not diffusion.")
