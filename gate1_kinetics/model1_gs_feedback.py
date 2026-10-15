"""
Model 1 — GS Feedback Inhibition: Wild-Type vs H195Y
Justifies selection of the H195Y mutant by quantifying its resistance to
product (glutamine) competitive inhibition at the glutamate active site.

All concentrations in mM. Time in minutes.
"""

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# --- Enzyme parameters (Fisher & Wray 2006, Tables 3–4) ---

ENZYMES = {
    "Wild type": dict(
        km_glu_mM=27.0, km_atp_mM=2.4, a_umol_min_mg=2.7, ic50_mM=2.4,
    ),
    "H195Y": dict(
        km_glu_mM=30.0, km_atp_mM=7.0, a_umol_min_mg=4.7, ic50_mM=115.0,
    ),
}

KM_NH3_MM = 0.4          # Nocek et al. (2008)
GLU_ASSAY_MM = 150.0     # Fisher & Wray (2006), Materials & Methods — used only in Cheng–Prusoff
GLU_FIXED_MM = 150.0     # operating glutamate
ATP_FIXED_MM = 3.5       # PPK2-buffered
E_MG_ML = 0.10           # comparison loading (same for both; ranking is loading-independent)
NH3_0_MM = 5.0
GLN_0_MM = 0.0
T_MAX_MIN = 60.0


def cheng_prusoff(ic50_mM, km_glu_mM, glu_assay_mM):
    """IC50 -> Ki, removing assay glutamate background (Cheng & Prusoff 1973)."""
    return ic50_mM / (1.0 + glu_assay_mM / km_glu_mM)


def gs_rate(nh3_mM, gln_mM, km_glu_mM, km_atp_mM, ki_mM, vmax_mM_min):
    """Three-substrate multiplicative rate with competitive glutamine inhibition on glutamate."""
    km_glu_app = km_glu_mM * (1.0 + gln_mM / ki_mM)
    f_atp = ATP_FIXED_MM / (km_atp_mM + ATP_FIXED_MM)
    f_nh3 = nh3_mM / (KM_NH3_MM + nh3_mM)
    f_glu = GLU_FIXED_MM / (km_glu_app + GLU_FIXED_MM)
    return vmax_mM_min * f_atp * f_nh3 * f_glu


def ode_system(t, y, km_glu_mM, km_atp_mM, ki_mM, vmax_mM_min):
    nh3, gln = y
    if nh3 < 0:
        nh3 = 0.0
    v = gs_rate(nh3, gln, km_glu_mM, km_atp_mM, ki_mM, vmax_mM_min)
    return [-v, +v]


# --- Compute Ki and print ---

results = {}
for name, p in ENZYMES.items():
    ki = cheng_prusoff(p["ic50_mM"], p["km_glu_mM"], GLU_ASSAY_MM)
    vmax = p["a_umol_min_mg"] * E_MG_ML  # µmol/min/mg × mg/mL = µmol/min/mL = mM/min
    results[name] = dict(ki=ki, vmax=vmax, **p)
    print(f"{name}: Ki = {ki:.3f} mM, Vmax = {vmax:.4f} mM/min")

ki_ratio = results["H195Y"]["ki"] / results["Wild type"]["ki"]
print(f"Ki ratio H195Y/WT = {ki_ratio:.1f}-fold")

# --- Rate vs accumulated glutamine (Figure 1) ---

gln_sweep = np.linspace(0, 50, 500)

fig1, ax1 = plt.subplots(figsize=(8, 5), constrained_layout=True)
colors = {"Wild type": "#E69F00", "H195Y": "#0072B2"}

for name, r in results.items():
    rates = [gs_rate(NH3_0_MM, g, r["km_glu_mM"], r["km_atp_mM"], r["ki"], r["vmax"])
             for g in gln_sweep]
    ax1.plot(gln_sweep, rates, label=name, color=colors[name], linewidth=2)

    # rate at 20 mM glutamine
    rate_20 = gs_rate(NH3_0_MM, 20.0, r["km_glu_mM"], r["km_atp_mM"], r["ki"], r["vmax"])
    rate_0 = gs_rate(NH3_0_MM, 0.0, r["km_glu_mM"], r["km_atp_mM"], r["ki"], r["vmax"])
    frac = rate_20 / rate_0 * 100
    print(f"{name} at [Gln]=20 mM: rate = {rate_20:.4f} mM/min ({frac:.0f}% of initial)")

ax1.set_xlabel("Accumulated glutamine (mM)", fontsize=12)
ax1.set_ylabel(r"$V_{\mathrm{GS}}$ (mM/min)", fontsize=12)
ax1.set_title("GS rate vs product accumulation\n"
              r"(5 mM NH$_3$, 150 mM glutamate, 3.5 mM ATP)", fontsize=12)
ax1.legend(frameon=True, fontsize=11)
ax1.set_xlim(0, 50)
ax1.set_ylim(bottom=0)
fig1.savefig("model1_fig1_rate_vs_glutamine.png", dpi=200, bbox_inches="tight")
plt.close(fig1)

# --- ODE solve: ammonia clearance over 60 min (Figure 2) ---

t_span = (0, T_MAX_MIN)
t_eval = np.linspace(0, T_MAX_MIN, 1000)

fig2, ax2 = plt.subplots(figsize=(8, 5), constrained_layout=True)

for name, r in results.items():
    sol = solve_ivp(
        ode_system, t_span, [NH3_0_MM, GLN_0_MM],
        args=(r["km_glu_mM"], r["km_atp_mM"], r["ki"], r["vmax"]),
        method="LSODA", t_eval=t_eval, rtol=1e-8, atol=1e-10,
    )
    nh3_final = sol.y[0, -1]
    gln_final = sol.y[1, -1]
    km_app_final = r["km_glu_mM"] * (1.0 + gln_final / r["ki"])
    clearance_pct = (1.0 - nh3_final / NH3_0_MM) * 100

    print(f"{name} at t=60 min: NH3 = {nh3_final:.3f} mM, cleared = {clearance_pct:.1f}%, "
          f"Km_app = {km_app_final:.1f} mM")

    ax2.plot(sol.t, sol.y[0], label=name, color=colors[name], linewidth=2)

ax2.set_xlabel("Time (min)", fontsize=12)
ax2.set_ylabel(r"[NH$_3$] (mM)", fontsize=12)
ax2.set_title(r"Ammonia clearance by wild type and H195Y GS", fontsize=12)
ax2.legend(frameon=True, fontsize=11)
ax2.set_xlim(0, T_MAX_MIN)
ax2.set_ylim(bottom=0)
fig2.savefig("model1_fig2_ammonia_clearance.png", dpi=200, bbox_inches="tight")
plt.close(fig2)

print("\nConclusion: H195Y is selected. Its main advantage is 52-fold weaker product "
      "binding (Ki = 19.2 vs 0.37 mM), not higher Vmax (only 1.7x).")
