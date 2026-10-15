"""
Model 7 — Batch Scrubbing Time and Farm Sizing
Computes hours and unit counts to bring a poultry shed from C0 to 20 ppm
using Ø 400 mm cartridges at 0.4 m/s, by solving:
    V_shed × dC/dt = -n × R(C)

R(C) is built from Model 6's plug-flow removal function, tabulated once
and interpolated for speed.

Concentrations in mol/m³ (gas phase). Time in seconds internally,
printed in hours.
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
from scipy.optimize import brentq
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

# --- Kinetic parameters ---

A_GS = 0.0783333
M_LOAD = 2.0
ETA_BEAD = 1.00
ETA_COAT = 0.985
KM_NH3_MOLM3 = 0.4
KM_GLU_MM = 30.0
KM_ATP_MM = 7.0
GLU_FIXED_MM = 150.0
ATP_FIXED_MM = 3.5

f_Glu = GLU_FIXED_MM / (KM_GLU_MM + GLU_FIXED_MM)
f_ATP = ATP_FIXED_MM / (KM_ATP_MM + ATP_FIXED_MM)


def vmax_eff_val(rho):
    return ETA_BEAD * ETA_COAT * rho * A_GS * M_LOAD * f_Glu * f_ATP


def wall_flux(C_gas, v_max_eff):
    jm = v_max_eff * L_COAT
    k_g = SH * D_NH3_AIR / D_h
    a_c = k_g * H_CC
    b_c = -(k_g * H_CC * C_gas - k_g * KM_NH3_MOLM3 - jm * H_CC)
    c_c = -k_g * KM_NH3_MOLM3 * C_gas
    disc = b_c**2 - 4 * a_c * c_c
    C_s = (-b_c + np.sqrt(max(disc, 0))) / (2 * a_c)
    c_liq = H_CC * C_s
    f_nh3 = c_liq / (KM_NH3_MOLM3 + c_liq)
    return v_max_eff * L_COAT * f_nh3


def removal_rate_one_unit(C_gas_molm3, D_mm, u_face, rho):
    """R in mol/s for one unit at given gas-phase concentration."""
    v_max_eff = vmax_eff_val(rho)
    D_m = D_mm * 1e-3
    N_cells = int(np.floor(PACK_FRACTION * np.pi * D_m**2 / 4 / A_cell))
    Q = u_face * np.pi * D_m**2 / 4
    u_ch = u_face / (PACK_FRACTION * eps_open)

    def dC_dx(x, C):
        if C[0] <= 0:
            return [0.0]
        j = wall_flux(C[0], v_max_eff)
        return [-(4.0 / (D_h * u_ch)) * j]

    sol = solve_ivp(dC_dx, [0, Z_STACK], [C_gas_molm3], method="RK45",
                    rtol=1e-10, atol=1e-14, max_step=Z_STACK / 50)
    C_out = max(sol.y[0, -1], 0.0)
    return Q * (C_gas_molm3 - C_out)


# --- Build interpolation table for R(C) ---

def build_R_interp(D_mm, u_face, rho):
    c_grid_ppm = np.linspace(1, 150, 300)
    c_grid_molm3 = c_grid_ppm * PPM_TO_MOLM3
    R_grid = np.array([removal_rate_one_unit(c, D_mm, u_face, rho) for c in c_grid_molm3])
    return interp1d(c_grid_molm3, R_grid, kind="cubic", fill_value="extrapolate")


print("Building R(C) interpolation tables...")
R_interp_nom = build_R_interp(400, 0.4, 0.5)
R_interp_lo = build_R_interp(400, 0.4, 0.3)
R_interp_hi = build_R_interp(400, 0.4, 0.7)


def batch_time_s(V_shed, C0_ppm, Cf_ppm, n_units, R_interp):
    """Solve V_shed dC/dt = -n R(C) until C = Cf. Returns time in seconds."""
    C0 = C0_ppm * PPM_TO_MOLM3
    Cf = Cf_ppm * PPM_TO_MOLM3

    def ode(t, C):
        if C[0] <= Cf:
            return [0.0]
        return [-n_units * R_interp(C[0]) / V_shed]

    event = lambda t, C: C[0] - Cf
    event.terminal = True
    event.direction = -1

    sol = solve_ivp(ode, [0, 200 * 3600], [C0], method="LSODA",
                    events=event, rtol=1e-8, atol=1e-12, max_step=600)

    if sol.t_events[0].size > 0:
        return sol.t_events[0][0]
    return sol.t[-1]


def batch_time_h(V_shed, C0_ppm, Cf_ppm, n_units, R_interp):
    return batch_time_s(V_shed, C0_ppm, Cf_ppm, n_units, R_interp) / 3600


# ============================================================
# Farm scenarios
# ============================================================

C_TARGET = 20.0

farms = [
    ("500 birds",   150, 60,  1),
    ("1,500 birds", 375, 80,  3),
    ("3,000 birds", 750, 100, 6),
]

print("\n=== Batch scrubbing times (Ø 400 mm, 0.4 m/s, ρ = 0.5) ===")
print(f"{'Farm':>20}  {'Shed':>6}  {'C0':>4}  {'Units':>5}  {'Time (h)':>8}  {'ρ=0.7':>6}  {'ρ=0.3':>6}")

farm_results = []
for label, V, C0, n in farms:
    t_nom = batch_time_h(V, C0, C_TARGET, n, R_interp_nom)
    t_hi = batch_time_h(V, C0, C_TARGET, n, R_interp_hi)
    t_lo = batch_time_h(V, C0, C_TARGET, n, R_interp_lo)

    # zero-order sanity check
    C_bar = (C0 + C_TARGET) / 2 * PPM_TO_MOLM3
    R_bar = n * R_interp_nom(C_bar)
    t_approx = V * (C0 - C_TARGET) * PPM_TO_MOLM3 / R_bar / 3600

    farm_results.append(dict(label=label, V=V, C0=C0, n=n,
                             t_nom=t_nom, t_hi=t_hi, t_lo=t_lo))
    print(f"{label:>20}  {V:>6}  {C0:>4}  {n:>5}  {t_nom:>8.1f}  {t_hi:>6.1f}  {t_lo:>6.1f}  "
          f"(approx {t_approx:.1f} h)")


# ============================================================
# Figure 1: concentration vs time for recommended deployments
# ============================================================

fig1, ax1 = plt.subplots(figsize=(9, 5.5), constrained_layout=True)
colors_farm = ["#0072B2", "#E69F00", "#D55E00"]

for i, (label, V, C0, n) in enumerate(farms):
    t_pts = np.linspace(0, 15 * 3600, 500)  # up to 15 h

    def ode_trace(t, C):
        if C[0] <= C_TARGET * PPM_TO_MOLM3:
            return [0.0]
        return [-n * R_interp_nom(C[0]) / V]

    sol = solve_ivp(ode_trace, [0, t_pts[-1]], [C0 * PPM_TO_MOLM3],
                    method="LSODA", t_eval=t_pts, rtol=1e-8, atol=1e-12)

    c_ppm = sol.y[0] / PPM_TO_MOLM3
    t_h = sol.t / 3600

    # ρ bands
    def ode_rho(t, C, R_int):
        if C[0] <= C_TARGET * PPM_TO_MOLM3:
            return [0.0]
        return [-n * R_int(C[0]) / V]

    sol_hi = solve_ivp(ode_rho, [0, t_pts[-1]], [C0 * PPM_TO_MOLM3],
                       args=(R_interp_hi,), method="LSODA", t_eval=t_pts, rtol=1e-8, atol=1e-12)
    sol_lo = solve_ivp(ode_rho, [0, t_pts[-1]], [C0 * PPM_TO_MOLM3],
                       args=(R_interp_lo,), method="LSODA", t_eval=t_pts, rtol=1e-8, atol=1e-12)

    ax1.fill_between(t_h, sol_hi.y[0] / PPM_TO_MOLM3, sol_lo.y[0] / PPM_TO_MOLM3,
                     color=colors_farm[i], alpha=0.12)
    ax1.plot(t_h, c_ppm, color=colors_farm[i], linewidth=2,
             label=f"{label} ({V} m³), {n} unit{'s' if n > 1 else ''}")

ax1.axhline(C_TARGET, color="gray", linestyle="--", linewidth=1, alpha=0.6)
ax1.text(14.5, C_TARGET + 2, "20 ppm target", fontsize=8, color="gray", ha="right")

ax1.set_xlabel("Time (h)", fontsize=11)
ax1.set_ylabel(r"Shed NH$_3$ (ppm)", fontsize=11)
ax1.set_title("Batch scrubbing with recommended units\n"
              r"(Ø 400 mm, 0.4 m/s, shaded: $\rho$ = 0.3–0.7)", fontsize=12)
ax1.legend(frameon=True, fontsize=9, loc="upper right")
ax1.set_xlim(0, 15)
ax1.set_ylim(0, 110)

fig1.savefig("model7_fig1_batch_scrubbing.png", dpi=200, bbox_inches="tight")
plt.close(fig1)


# ============================================================
# Figure 2: batch time vs starting concentration (1,500-bird)
# ============================================================

V_1500 = 375
c0_sweep = np.linspace(25, 150, 200)
unit_options = [2, 3, 4, 5]
unit_colors = ["#56B4E9", "#0072B2", "#E69F00", "#D55E00"]

fig2, ax2 = plt.subplots(figsize=(9, 5.5), constrained_layout=True)

crossings = {}
for n_u, col in zip(unit_options, unit_colors):
    times = [batch_time_h(V_1500, c0, C_TARGET, n_u, R_interp_nom) for c0 in c0_sweep]
    ax2.plot(c0_sweep, times, color=col, linewidth=2, label=f"{n_u} units")

    # find crossing with 8 h line
    times_arr = np.array(times)
    above_8 = c0_sweep[times_arr > 8.0]
    if len(above_8) > 0 and len(above_8) < len(c0_sweep):
        try:
            t_interp_func = interp1d(c0_sweep, times_arr - 8.0, kind="cubic")
            c_cross = brentq(t_interp_func, c0_sweep[0], c0_sweep[-1])
            crossings[n_u] = c_cross
            ax2.plot(c_cross, 8.0, 'o', color=col, markersize=7, zorder=5)
            ax2.annotate(f"{c_cross:.0f} ppm", xy=(c_cross, 8.0),
                         xytext=(c_cross + 8, 8.8),
                         arrowprops=dict(arrowstyle="->", color=col),
                         fontsize=8, color=col)
        except ValueError:
            pass

ax2.axhline(8.0, color="gray", linestyle="--", linewidth=1, alpha=0.6)
ax2.text(148, 8.3, "8-hour run limit", fontsize=8, color="gray", ha="right")

ax2.set_xlabel(r"Starting NH$_3$ (ppm)", fontsize=11)
ax2.set_ylabel("Batch time to 20 ppm (h)", fontsize=11)
ax2.set_title("Batch time vs starting concentration\n"
              "(1,500-bird farm, 375 m³, Ø 400 mm units)", fontsize=12)
ax2.legend(frameon=True, fontsize=10, loc="upper left")
ax2.set_xlim(20, 150)
ax2.set_ylim(0, 20)

fig2.savefig("model7_fig2_time_vs_C0.png", dpi=200, bbox_inches="tight")
plt.close(fig2)

print("\n8-hour crossing concentrations (1,500-bird farm):")
for n_u, c_cross in crossings.items():
    print(f"  {n_u} units: up to {c_cross:.0f} ppm")


# ============================================================
# Figure 3: Sensitivity to ATP and glutamate (1,500-bird, 3 units)
# ============================================================

def build_R_interp_chem(rho, glu_mm, atp_mm):
    """Build R(C) table with modified f_Glu and f_ATP."""
    f_g = glu_mm / (KM_GLU_MM + glu_mm)
    f_a = atp_mm / (KM_ATP_MM + atp_mm)
    v = ETA_BEAD * ETA_COAT * rho * A_GS * M_LOAD * f_g * f_a

    c_grid = np.linspace(1, 150, 200) * PPM_TO_MOLM3
    R_grid = []
    for c in c_grid:
        u_ch = 0.4 / (PACK_FRACTION * eps_open)
        D_m = 0.4
        N_cells = int(np.floor(PACK_FRACTION * np.pi * D_m**2 / 4 / A_cell))
        Q = 0.4 * np.pi * D_m**2 / 4

        def dC_dx(x, C):
            if C[0] <= 0:
                return [0.0]
            jm = v * L_COAT
            k_g = SH * D_NH3_AIR / D_h
            a_co = k_g * H_CC
            b_co = -(k_g * H_CC * C[0] - k_g * KM_NH3_MOLM3 - jm * H_CC)
            c_co = -k_g * KM_NH3_MOLM3 * C[0]
            disc = b_co**2 - 4 * a_co * c_co
            C_s = (-b_co + np.sqrt(max(disc, 0))) / (2 * a_co)
            c_liq = H_CC * C_s
            f_nh3 = c_liq / (KM_NH3_MOLM3 + c_liq)
            j_val = v * L_COAT * f_nh3
            return [-(4.0 / (D_h * u_ch)) * j_val]

        sol = solve_ivp(dC_dx, [0, Z_STACK], [c], method="RK45",
                        rtol=1e-10, atol=1e-14, max_step=Z_STACK / 50)
        C_out = max(sol.y[0, -1], 0.0)
        R_grid.append(Q * (c - C_out))

    return interp1d(c_grid, R_grid, kind="cubic", fill_value="extrapolate")


chem_cases = [
    ("Baseline (3.5 mM ATP, 150 mM Glu)", 150.0, 3.5, "#0072B2"),
    ("7 mM ATP", 150.0, 7.0, "#009E73"),
    ("1.75 mM ATP", 150.0, 1.75, "#D55E00"),
    ("300 mM Glu", 300.0, 3.5, "#CC79A7"),
]

fig3, ax3 = plt.subplots(figsize=(9, 5.5), constrained_layout=True)

print("\nSensitivity (1,500-bird farm, 3 units, 80→20 ppm):")
for label, glu, atp, col in chem_cases:
    R_int = build_R_interp_chem(0.5, glu, atp)
    t = batch_time_h(V_1500, 80, C_TARGET, 3, R_int)
    print(f"  {label}: {t:.1f} h")

    # time trace
    t_pts = np.linspace(0, 15 * 3600, 300)

    def ode_chem(t_val, C, R_i=R_int):
        if C[0] <= C_TARGET * PPM_TO_MOLM3:
            return [0.0]
        return [-3 * R_i(C[0]) / V_1500]

    sol = solve_ivp(ode_chem, [0, t_pts[-1]], [80 * PPM_TO_MOLM3],
                    method="LSODA", t_eval=t_pts, rtol=1e-8, atol=1e-12)
    ax3.plot(sol.t / 3600, sol.y[0] / PPM_TO_MOLM3, color=col, linewidth=2, label=label)

ax3.axhline(C_TARGET, color="gray", linestyle="--", linewidth=1, alpha=0.5)
ax3.set_xlabel("Time (h)", fontsize=11)
ax3.set_ylabel(r"Shed NH$_3$ (ppm)", fontsize=11)
ax3.set_title("Sensitivity of batch scrubbing to ATP and glutamate\n"
              "(1,500-bird farm, 3 units)", fontsize=12)
ax3.legend(frameon=True, fontsize=9, loc="upper right")
ax3.set_xlim(0, 15)
ax3.set_ylim(0, 90)

fig3.savefig("model7_fig3_sensitivity.png", dpi=200, bbox_inches="tight")
plt.close(fig3)


# ============================================================
# Figure 4: Units needed within 10 h vs diameter
# ============================================================

diameters = [200, 250, 300, 400, 500]

farm_specs = [
    ("500-bird farm (150 m³)", 150, 60),
    ("1,500-bird farm (375 m³)", 375, 80),
    ("3,000-bird farm (750 m³)", 750, 100),
]
farm_colors = ["#0072B2", "#E69F00", "#D55E00"]

print("\nUnits needed to reach 20 ppm within 10 h:")
print(f"{'D (mm)':>8}", end="")
for lbl, _, _ in farm_specs:
    print(f"  {lbl:>25}", end="")
print()

units_needed = {lbl: [] for lbl, _, _ in farm_specs}
for D_mm in diameters:
    R_int_d = build_R_interp(D_mm, 0.4, 0.5)
    row = f"{D_mm:>8}"
    for lbl, V, C0 in farm_specs:
        for n in range(1, 100):
            t = batch_time_h(V, C0, C_TARGET, n, R_int_d)
            if t <= 10.0:
                units_needed[lbl].append(n)
                row += f"  {n:>25}"
                break
        else:
            units_needed[lbl].append(99)
            row += f"  {'99+':>25}"
    print(row)

fig4, ax4 = plt.subplots(figsize=(9, 5.5), constrained_layout=True)
x_pos = np.arange(len(diameters))
width = 0.25

for i, (lbl, _, _) in enumerate(farm_specs):
    ax4.bar(x_pos + i * width - width, units_needed[lbl], width,
            color=farm_colors[i], label=lbl)
    for j, n in enumerate(units_needed[lbl]):
        ax4.text(x_pos[j] + i * width - width, n + 0.3, str(n),
                 ha="center", va="bottom", fontsize=8, color=farm_colors[i])

ax4.set_xlabel("Cartridge diameter (mm)", fontsize=11)
ax4.set_ylabel("Units needed", fontsize=11)
ax4.set_title("Units needed to reach 20 ppm within 10 h\n(Ø 400 mm, 0.4 m/s, ρ = 0.5)", fontsize=12)
ax4.set_xticks(x_pos)
ax4.set_xticklabels([str(d) for d in diameters])
ax4.legend(frameon=True, fontsize=9, loc="upper right")
ax4.set_ylim(0, max(max(v) for v in units_needed.values()) * 1.3)

fig4.savefig("model7_fig4_units_vs_diameter.png", dpi=200, bbox_inches="tight")
plt.close(fig4)
