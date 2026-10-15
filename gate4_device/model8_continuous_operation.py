"""
Model 8 — Continuous Operation Check
Tests whether running cartridges continuously to hold 20 ppm is
practical, and compares against batch mode (Model 7).

N_units = (S - Q_vent × C_target) / R(C_target)
where S = Q_1 × C_0 (emission at 1 ACH).
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

# --- Physical constants ---

R_GAS = 8.314462
T_K = 298.15
P_ATM = 101325.0
PPM_TO_MOLM3 = P_ATM / (R_GAS * T_K) * 1e-6

H_CC = 1443.0
D_NH3_AIR = 2.3e-5
SH = 3.0
MU_AIR = 1.84e-5

# --- Scaffold geometry ---

L_COAT = 24e-6
R_CELL = 1.5e-3
T_WALL = 0.8e-3
PACK_FRACTION = 0.92
Z_STACK = 14 * 0.020

r_open = R_CELL - L_COAT
A_open = 2 * np.sqrt(3) * r_open**2
A_cell = 2 * np.sqrt(3) * (R_CELL + T_WALL / 2)**2
P_wet = 4 * np.sqrt(3) * r_open
D_h = 2 * r_open
eps_open = A_open / A_cell

# --- Kinetics ---

A_GS = 0.0783333
M_LOAD = 2.0
ETA_BEAD = 1.00
ETA_COAT = 0.985
RHO_RETAINED = 0.5
KM_NH3_MOLM3 = 0.4
KM_GLU_MM = 30.0
KM_ATP_MM = 7.0
GLU_FIXED_MM = 150.0
ATP_FIXED_MM = 3.5

f_Glu = GLU_FIXED_MM / (KM_GLU_MM + GLU_FIXED_MM)
f_ATP = ATP_FIXED_MM / (KM_ATP_MM + ATP_FIXED_MM)
V_MAX_EFF = ETA_BEAD * ETA_COAT * RHO_RETAINED * A_GS * M_LOAD * f_Glu * f_ATP


def removal_rate_one_unit(C_gas_molm3):
    """R in mol/s for one Ø 400 mm unit at 0.4 m/s."""
    u_face = 0.4
    D_m = 0.4
    N_cells = int(np.floor(PACK_FRACTION * np.pi * D_m**2 / 4 / A_cell))
    Q = u_face * np.pi * D_m**2 / 4
    u_ch = u_face / (PACK_FRACTION * eps_open)

    def dC_dx(x, C):
        if C[0] <= 0:
            return [0.0]
        jm = V_MAX_EFF * L_COAT
        k_g = SH * D_NH3_AIR / D_h
        a_c = k_g * H_CC
        b_c = -(k_g * H_CC * C[0] - k_g * KM_NH3_MOLM3 - jm * H_CC)
        c_c = -k_g * KM_NH3_MOLM3 * C[0]
        disc = b_c**2 - 4 * a_c * c_c
        C_s = (-b_c + np.sqrt(max(disc, 0))) / (2 * a_c)
        c_liq = H_CC * C_s
        f_nh3 = c_liq / (KM_NH3_MOLM3 + c_liq)
        return [-(4.0 / (D_h * u_ch)) * V_MAX_EFF * L_COAT * f_nh3]

    sol = solve_ivp(dC_dx, [0, Z_STACK], [C_gas_molm3], method="RK45",
                    rtol=1e-10, atol=1e-14, max_step=Z_STACK / 50)
    C_out = max(sol.y[0, -1], 0.0)
    return Q * (C_gas_molm3 - C_out)


# R at 20 ppm (the target)
C_TARGET_PPM = 20.0
C_TARGET_MOLM3 = C_TARGET_PPM * PPM_TO_MOLM3
R_at_target = removal_rate_one_unit(C_TARGET_MOLM3)
R_at_target_umol_min = R_at_target * 1e6 * 60

print(f"R(20 ppm) = {R_at_target_umol_min:.0f} µmol/min per unit")

# --- Farm scenarios ---

farms = [
    ("500-bird farm (150 m³)",   150.0,  60.0),
    ("1,500-bird farm (375 m³)", 375.0,  80.0),
    ("3,000-bird farm (750 m³)", 750.0, 100.0),
]

print(f"\n{'Farm':>30}  {'Emission (g/day)':>16}  {'Units at 1 ACH':>15}  {'Zero-crossing ACH':>18}")

farm_results = []
for label, V_shed, C0_ppm in farms:
    Q_1 = V_shed / 3600.0  # m³/s at 1 ACH
    S = Q_1 * C0_ppm * PPM_TO_MOLM3  # mol/s emission

    # emission in g/day
    emission_g_day = S * 17.031 * 86400  # NH3 MW = 17.031 g/mol

    # units at 1 ACH
    N_units_1ach = (S - Q_1 * C_TARGET_MOLM3) / R_at_target
    N_units_1ach_ceil = int(np.ceil(N_units_1ach))

    # ACH at which ventilation alone holds 20 ppm: ACH = C0 / C_target
    ach_zero = C0_ppm / C_TARGET_PPM

    farm_results.append(dict(label=label, V=V_shed, C0=C0_ppm,
                              emission_g=emission_g_day, N_1ach=N_units_1ach_ceil,
                              ach_zero=ach_zero))

    print(f"{label:>30}  {emission_g_day:>16.0f}  {N_units_1ach_ceil:>15d}  {ach_zero:>18.0f}")

# batch-mode comparison
batch_units = [1, 3, 6]
print(f"\nContinuous operation at 1 ACH requires {farm_results[0]['N_1ach']}–{farm_results[2]['N_1ach']} "
      f"units — versus {batch_units[0]}–{batch_units[2]} units in batch mode (Model 7).")
for fr, bu in zip(farm_results, batch_units):
    ratio = fr["N_1ach"] / bu
    print(f"  {fr['label']}: {fr['N_1ach']} vs {bu} units ({ratio:.0f}× reduction with batch)")


# ============================================================
# Figure: units needed vs ventilation rate (ACH)
# ============================================================

fig, ax = plt.subplots(figsize=(9, 5.5), constrained_layout=True)
colors = ["#0072B2", "#E69F00", "#D55E00"]

for i, fr in enumerate(farm_results):
    V_shed = fr["V"]
    C0 = fr["C0"]
    ach_zero = fr["ach_zero"]

    ach_sweep = np.linspace(1.0, ach_zero, 200)
    n_units_sweep = []
    for ach in ach_sweep:
        Q_vent = V_shed * ach / 3600.0  # m³/s
        S = (V_shed / 3600.0) * C0 * PPM_TO_MOLM3  # emission (based on 1 ACH)
        excess = S - Q_vent * C_TARGET_MOLM3
        if excess <= 0:
            n_units_sweep.append(0)
        else:
            n_units_sweep.append(excess / R_at_target)

    ax.plot(ach_sweep, n_units_sweep, color=colors[i], linewidth=2, label=fr["label"])

    # mark 1 ACH
    ax.plot(1.0, fr["N_1ach"], 'o', color=colors[i], markersize=7, zorder=5)
    ax.annotate(f"{fr['N_1ach']} units",
                xy=(1.0, fr["N_1ach"]),
                xytext=(1.3, fr["N_1ach"] + 3),
                arrowprops=dict(arrowstyle="->", color=colors[i]),
                fontsize=9, color=colors[i])

    # mark zero crossing
    ax.plot(ach_zero, 0, 'v', color=colors[i], markersize=7, zorder=5)
    ax.annotate(f"{ach_zero:.0f} ACH",
                xy=(ach_zero, 0),
                xytext=(ach_zero, 5),
                arrowprops=dict(arrowstyle="->", color=colors[i]),
                fontsize=9, color=colors[i], ha="center")

ax.set_xlabel("Ventilation rate (ACH)", fontsize=11)
ax.set_ylabel("Units needed to hold 20 ppm", fontsize=11)
ax.set_title("Continuous operation: units needed vs ventilation\n"
             "(Ø 400 mm, 0.4 m/s, ρ = 0.5)", fontsize=12)
ax.legend(frameon=True, fontsize=10, loc="upper right")
ax.set_xlim(0.8, 6)
ax.set_ylim(bottom=0)

fig.savefig("model8_continuous_operation.png", dpi=200, bbox_inches="tight")
plt.close(fig)

print("\nConclusion: Continuous operation at realistic ventilation needs 7–70 units — "
      "not affordable for the target farms. Batch scrubbing during peaks reaches 20 ppm "
      "with 1–6 units (a 7–12× reduction). Good ventilation is still the cheapest first step.")
