"""
Model 6 — Single-Pass Removal Rate of One Cartridge
Combines enzyme kinetics with plug-flow integration along the channel
to give the actual removal rate, respecting the constraint that removal
cannot exceed what the airflow delivers.

Length in m, time in s, concentration in mol/m³ (gas phase).
"""

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# --- Physical constants ---

R_GAS = 8.314462
T_K = 298.15
P_ATM = 101325.0
PPM_TO_MOLM3 = P_ATM / (R_GAS * T_K) * 1e-6   # 4.087e-5 mol/m³ per ppm

H_CC = 1443.0
D_NH3_AIR = 2.3e-5          # m²/s
SH = 3.0
MU_AIR = 1.84e-5            # Pa·s
RHO_AIR = 1.184              # kg/m³

# --- Scaffold geometry ---

L_COAT = 24e-6
R_CELL = 1.5e-3
T_WALL = 0.8e-3
PACK_FRACTION = 0.92
Z_TRAY = 0.020
N_TRAYS = 14
Z_STACK = N_TRAYS * Z_TRAY   # 0.280 m

r_open = R_CELL - L_COAT
r_out = R_CELL + T_WALL / 2
A_open = 2 * np.sqrt(3) * r_open**2
A_cell = 2 * np.sqrt(3) * r_out**2
P_wet = 4 * np.sqrt(3) * r_open
D_h = 2 * r_open
eps_open = A_open / A_cell

# --- Kinetic parameters ---

A_GS = 0.0783333             # mol s⁻¹ kg⁻¹
M_LOAD = 2.0                 # kg/m³
ETA_BEAD = 1.00
ETA_COAT = 0.985
RHO_RETAINED = 0.5

KM_GLU_MM = 30.0
KM_ATP_MM = 7.0
KM_NH3_MM = 0.4
GLU_FIXED_MM = 150.0
ATP_FIXED_MM = 3.5

f_Glu = GLU_FIXED_MM / (KM_GLU_MM + GLU_FIXED_MM)   # 0.833
f_ATP = ATP_FIXED_MM / (KM_ATP_MM + ATP_FIXED_MM)    # 0.333

KM_NH3_MOLM3 = KM_NH3_MM    # mM = mol/m³

def vmax_eff(rho):
    return ETA_BEAD * ETA_COAT * rho * A_GS * M_LOAD * f_Glu * f_ATP

V_MAX_EFF = vmax_eff(RHO_RETAINED)     # mol/(s·m³) in coating volume
V_MAX_EFF_LO = vmax_eff(0.3)
V_MAX_EFF_HI = vmax_eff(0.7)

j_m_nom = V_MAX_EFF * L_COAT           # mol/(m²·s), saturated wall flux
k_g = SH * D_NH3_AIR / D_h             # m/s, gas-film coefficient


def wall_flux(C_gas, v_max_eff_local):
    """Solve for C_s from the gas-film + reaction balance, then compute j.
    k_g (C - C_s) = v_max_eff * L_coat * H_cc*C_s / (Km_NH3 + H_cc*C_s)
    Rearranged as a quadratic in C_s (positive root).
    """
    jm = v_max_eff_local * L_COAT
    # coefficients: k_g*H_cc * Cs^2 - (k_g*H_cc*C - k_g*Km_NH3 - jm*H_cc)*Cs - k_g*Km_NH3*C = 0
    a_coeff = k_g * H_CC
    b_coeff = -(k_g * H_CC * C_gas - k_g * KM_NH3_MOLM3 - jm * H_CC)
    c_coeff = -k_g * KM_NH3_MOLM3 * C_gas

    disc = b_coeff**2 - 4 * a_coeff * c_coeff
    C_s = (-b_coeff + np.sqrt(np.maximum(disc, 0))) / (2 * a_coeff)

    c_liq = H_CC * C_s
    f_nh3 = c_liq / (KM_NH3_MOLM3 + c_liq)
    j = v_max_eff_local * L_COAT * f_nh3
    return j, C_s


def plug_flow_outlet(C_in, u_face, D_cart, v_max_eff_local):
    """Integrate dC/dx = -(4/(D_h·u_ch)) · j(C) along x ∈ [0, z_stack]."""
    u_ch = u_face / (PACK_FRACTION * eps_open)

    def dC_dx(x, C):
        if C[0] <= 0:
            return [0.0]
        j, _ = wall_flux(C[0], v_max_eff_local)
        return [-(4.0 / (D_h * u_ch)) * j]

    sol = solve_ivp(dC_dx, [0, Z_STACK], [C_in], method="RK45",
                    rtol=1e-10, atol=1e-14, max_step=Z_STACK / 100)
    return max(sol.y[0, -1], 0.0), u_ch


def cartridge_metrics(C_in_ppm, u_face, D_cart_mm, v_max_eff_local):
    """Full set of metrics for one cartridge at given conditions."""
    C_in = C_in_ppm * PPM_TO_MOLM3
    D_cart = D_cart_mm * 1e-3

    N_cells = int(np.floor(PACK_FRACTION * np.pi * D_cart**2 / 4 / A_cell))
    Q = u_face * np.pi * D_cart**2 / 4  # m³/s
    TCWA = N_cells * P_wet * Z_STACK     # m²

    C_out, u_ch = plug_flow_outlet(C_in, u_face, D_cart, v_max_eff_local)

    R_actual = Q * (C_in - C_out)        # mol/s
    R_del = Q * C_in                     # mol/s, delivery limit
    j_in, _ = wall_flux(C_in, v_max_eff_local)
    R_kin = j_in * TCWA                  # mol/s, kinetic ceiling

    X = 1.0 - C_out / C_in if C_in > 0 else 0.0

    Re = RHO_AIR * u_ch * D_h / MU_AIR
    dP = 32 * MU_AIR * u_ch * Z_STACK / D_h**2

    # film flux loss
    j_no_film, _ = wall_flux(C_in, v_max_eff_local)
    c_liq_direct = H_CC * C_in
    f_nh3_direct = c_liq_direct / (KM_NH3_MOLM3 + c_liq_direct)
    j_direct = v_max_eff_local * L_COAT * f_nh3_direct
    film_loss_pct = (1.0 - j_no_film / j_direct) * 100 if j_direct > 0 else 0.0

    return dict(
        N_cells=N_cells, Q_m3s=Q, Q_m3h=Q * 3600, TCWA=TCWA,
        C_out_ppm=C_out / PPM_TO_MOLM3,
        R_umol_min=R_actual * 1e6 * 60, R_del_umol_min=R_del * 1e6 * 60,
        R_kin_umol_min=R_kin * 1e6 * 60,
        X=X, u_ch=u_ch, Re=Re, dP=dP, film_loss_pct=film_loss_pct,
    )


# ============================================================
# Print expected results: D = 400 mm, u_face = 0.4 m/s, 50 ppm
# ============================================================

m = cartridge_metrics(50.0, 0.4, 400, V_MAX_EFF)
print("=== D = 400 mm, u_face = 0.4 m/s, 50 ppm ===")
print(f"N_cells       = {m['N_cells']}")
print(f"TCWA          = {m['TCWA']:.2f} m²")
print(f"Q             = {m['Q_m3h']:.1f} m³/h")
print(f"u_ch          = {m['u_ch']:.2f} m/s")
print(f"Re            = {m['Re']:.0f} (laminar)")
print(f"k_g           = {k_g:.4f} m/s")
print(f"R             = {m['R_umol_min']:.0f} µmol/min")
print(f"R_kin         = {m['R_kin_umol_min']:.0f} µmol/min")
print(f"R_del         = {m['R_del_umol_min']:.0f} µmol/min  ({m['R_del_umol_min']/m['R_umol_min']:.1f}× removal)")
print(f"Single-pass X = {m['X']*100:.1f}%")
print(f"C_out         = {m['C_out_ppm']:.1f} ppm")
print(f"Film loss     = {m['film_loss_pct']:.2f}%")
print(f"ΔP scaffold   = {m['dP']:.1f} Pa")
print("The unit is kinetics-limited, not supply-limited.")

# pipe fit note
print(f"\nNote: a nominal 16″ PVC pipe has OD ~406 mm; inner diameter is below 400 mm. "
      f"A 385 mm disc would have ~{(1 - (385/400)**2)*100:.0f}% fewer cells.")

# pressure budget
print(f"\nPressure budget: scaffold {m['dP']:.1f} Pa + pre-filter 20–50 Pa "
      f"+ humidifier 30–100 Pa = {m['dP']+20:.0f}–{m['dP']+150:.0f} Pa total")

# R_del = R_kin crossover velocity
u_sweep_cross = np.linspace(0.01, 0.15, 1000)
for u in u_sweep_cross:
    mc = cartridge_metrics(50.0, u, 400, V_MAX_EFF)
    if mc["R_del_umol_min"] >= mc["R_kin_umol_min"]:
        print(f"Delivery = kinetic ceiling at u_face ≈ {u:.3f} m/s")
        break

# within 0.3–0.5 range, check closeness to kinetic ceiling
for u_test in [0.3, 0.5]:
    mc = cartridge_metrics(50.0, u_test, 400, V_MAX_EFF)
    gap = (1.0 - mc["R_umol_min"] / mc["R_kin_umol_min"]) * 100
    print(f"At u_face = {u_test}: R = {mc['R_umol_min']:.0f}, gap to R_kin = {gap:.1f}%")

# --- Concentration band ---

print("\nConcentration band (ρ = 0.5 nominal, 0.3/0.7 band):")
for c_ppm in [20, 50, 100]:
    c_liq = H_CC * c_ppm * PPM_TO_MOLM3
    f_nh3 = c_liq / (KM_NH3_MOLM3 + c_liq)
    r_nom = cartridge_metrics(c_ppm, 0.4, 400, V_MAX_EFF)["R_umol_min"]
    r_lo = cartridge_metrics(c_ppm, 0.4, 400, V_MAX_EFF_LO)["R_umol_min"]
    r_hi = cartridge_metrics(c_ppm, 0.4, 400, V_MAX_EFF_HI)["R_umol_min"]
    print(f"  {c_ppm:3d} ppm: R = {r_nom:.0f} ({r_lo:.0f}–{r_hi:.0f}) µmol/min, f_NH3 = {f_nh3*100:.1f}%")

# --- Diameter table ---

print("\nDiameter table (u_face = 0.4 m/s, ρ = 0.5):")
print(f"{'D (mm)':>8} {'Cells':>7} {'TCWA (m²)':>10} {'R@50ppm':>10} {'Q (m³/h)':>10} {'PETG (kg)':>10}")
for D_mm in [200, 250, 300, 400, 500]:
    mc = cartridge_metrics(50.0, 0.4, D_mm, V_MAX_EFF)
    # PETG mass: approximate as wall volume × density (1.27 g/cm³)
    D_m = D_mm * 1e-3
    n_cells = mc["N_cells"]
    wall_vol = n_cells * (A_cell - A_open) * Z_STACK  # m³
    petg_kg = wall_vol * 1270  # kg
    print(f"{D_mm:8d} {n_cells:7d} {mc['TCWA']:10.2f} {mc['R_umol_min']:10.0f} "
          f"{mc['Q_m3h']:10.1f} {petg_kg:10.1f}")


# ============================================================
# Figures
# ============================================================

# --- Figure 1: (a) Removal rate and (b) single-pass capture vs face velocity ---

u_sweep = np.linspace(0.01, 0.8, 300)
R_arr = np.zeros_like(u_sweep)
R_del_arr = np.zeros_like(u_sweep)
R_kin_arr = np.zeros_like(u_sweep)
X_arr = np.zeros_like(u_sweep)

for i, u in enumerate(u_sweep):
    mc = cartridge_metrics(50.0, u, 400, V_MAX_EFF)
    R_arr[i] = mc["R_umol_min"]
    R_del_arr[i] = mc["R_del_umol_min"]
    R_kin_arr[i] = mc["R_kin_umol_min"]
    X_arr[i] = mc["X"] * 100

fig1, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(14, 5.5), constrained_layout=True)

ax_a.plot(u_sweep, R_arr, color="#0072B2", linewidth=2.5, label="Removal rate")
ax_a.plot(u_sweep, R_del_arr, color="#D55E00", linewidth=1.5, linestyle="--", label="Delivery (supply)")
ax_a.plot(u_sweep, R_kin_arr, color="#009E73", linewidth=1.5, linestyle=":", label="Kinetic ceiling")

# operating point
mc_op = cartridge_metrics(50.0, 0.4, 400, V_MAX_EFF)
ax_a.plot(0.4, mc_op["R_umol_min"], 'o', color="#0072B2", markersize=8, zorder=5)
ax_a.annotate(f"u = 0.4 m/s\nR = {mc_op['R_umol_min']:.0f} µmol/min",
              xy=(0.4, mc_op["R_umol_min"]),
              xytext=(0.5, mc_op["R_umol_min"] - 150),
              arrowprops=dict(arrowstyle="->", color="#0072B2"),
              fontsize=9, color="#0072B2")

ax_a.set_xlabel("Face velocity (m/s)", fontsize=11)
ax_a.set_ylabel("Rate (µmol/min)", fontsize=11)
ax_a.set_title("(a) Removal rate vs face velocity\n(Ø 400 mm, 50 ppm)", fontsize=12)
ax_a.legend(frameon=True, fontsize=9, loc="center right")
ax_a.set_xlim(0, 0.8)
ax_a.set_ylim(0, max(R_del_arr.max(), R_kin_arr[-1]) * 1.1)

ax_b.plot(u_sweep, X_arr, color="#0072B2", linewidth=2.5)
ax_b.plot(0.4, mc_op["X"] * 100, 'o', color="#0072B2", markersize=8, zorder=5)
ax_b.annotate(f"u = 0.4 m/s\nX = {mc_op['X']*100:.1f}%",
              xy=(0.4, mc_op["X"] * 100),
              xytext=(0.55, mc_op["X"] * 100 + 15),
              arrowprops=dict(arrowstyle="->", color="#0072B2"),
              fontsize=9, color="#0072B2")
ax_b.set_xlabel("Face velocity (m/s)", fontsize=11)
ax_b.set_ylabel("Single-pass capture (%)", fontsize=11)
ax_b.set_title("(b) Single-pass capture vs face velocity", fontsize=12)
ax_b.set_xlim(0, 0.8)
ax_b.set_ylim(0, 105)

fig1.savefig("model6_fig1_removal_vs_velocity.png", dpi=200, bbox_inches="tight")
plt.close(fig1)


# --- Figure 2: Removal rate vs shed concentration, with ρ band ---

c_sweep_ppm = np.linspace(5, 120, 200)
R_nom = np.array([cartridge_metrics(c, 0.4, 400, V_MAX_EFF)["R_umol_min"] for c in c_sweep_ppm])
R_lo = np.array([cartridge_metrics(c, 0.4, 400, V_MAX_EFF_LO)["R_umol_min"] for c in c_sweep_ppm])
R_hi = np.array([cartridge_metrics(c, 0.4, 400, V_MAX_EFF_HI)["R_umol_min"] for c in c_sweep_ppm])

fig2, ax2 = plt.subplots(figsize=(8, 5.5), constrained_layout=True)
ax2.fill_between(c_sweep_ppm, R_lo, R_hi, color="#0072B2", alpha=0.15, label=r"$\rho$ = 0.3–0.7")
ax2.plot(c_sweep_ppm, R_nom, color="#0072B2", linewidth=2.5, label=r"$\rho$ = 0.5")
ax2.set_xlabel(r"Shed NH$_3$ (ppm)", fontsize=11)
ax2.set_ylabel("Removal rate (µmol/min)", fontsize=11)
ax2.set_title(r"Removal rate of one Ø 400 mm unit at 0.4 m/s", fontsize=12)
ax2.legend(frameon=True, fontsize=10)
ax2.set_xlim(0, 120)
ax2.set_ylim(0, R_hi.max() * 1.1)

fig2.savefig("model6_fig2_removal_vs_concentration.png", dpi=200, bbox_inches="tight")
plt.close(fig2)


# --- Figure 3: Removal and airflow vs diameter ---

diameters = [200, 250, 300, 400, 500]
R_by_d = []
Q_by_d = []
for d in diameters:
    mc = cartridge_metrics(50.0, 0.4, d, V_MAX_EFF)
    R_by_d.append(mc["R_umol_min"])
    Q_by_d.append(mc["Q_m3h"])

fig3, ax3 = plt.subplots(figsize=(8, 5.5), constrained_layout=True)
x_pos = np.arange(len(diameters))
width = 0.35

ax3.bar(x_pos - width / 2, R_by_d, width, color="#0072B2", label="R at 50 ppm (µmol/min)")
ax3_twin = ax3.twinx()
ax3_twin.bar(x_pos + width / 2, Q_by_d, width, color="#E69F00", label="Q (m³/h)")

ax3.set_xlabel("Cartridge diameter (mm)", fontsize=11)
ax3.set_ylabel("Removal rate (µmol/min)", fontsize=11, color="#0072B2")
ax3_twin.set_ylabel("Airflow Q (m³/h)", fontsize=11, color="#E69F00")
ax3.set_title("Removal rate and airflow by cartridge diameter\n(u = 0.4 m/s, ρ = 0.5)", fontsize=12)
ax3.set_xticks(x_pos)
ax3.set_xticklabels([str(d) for d in diameters])
ax3.tick_params(axis="y", colors="#0072B2")
ax3_twin.tick_params(axis="y", colors="#E69F00")

lines_3 = ax3.containers + ax3_twin.containers
labels_3 = ["R at 50 ppm (µmol/min)", "Q (m³/h)"]
ax3.legend(lines_3, labels_3, frameon=True, fontsize=9, loc="upper left")

fig3.savefig("model6_fig3_diameter_comparison.png", dpi=200, bbox_inches="tight")
plt.close(fig3)


# --- Figure 4: Pressure drop budget ---

fig4, ax4 = plt.subplots(figsize=(8, 3), constrained_layout=True)
categories = ["Scaffold", "Pre-filter\n(estimate)", "Humidifier\n(estimate)"]
vals_lo = [m["dP"], 20, 30]
vals_hi = [m["dP"], 50, 100]

left = 0
for i, (cat, lo, hi) in enumerate(zip(categories, vals_lo, vals_hi)):
    color = ["#0072B2", "#E69F00", "#009E73"][i]
    bar_len = (lo + hi) / 2
    ax4.barh(0, bar_len, left=left, color=color, edgecolor="white", height=0.5, label=cat)
    ax4.text(left + bar_len / 2, 0, f"{lo:.0f}–{hi:.0f} Pa" if lo != hi else f"{lo:.0f} Pa",
             ha="center", va="center", fontsize=9, color="white", fontweight="bold")
    left += bar_len

ax4.set_xlabel("Pressure drop (Pa)", fontsize=11)
ax4.set_title("Pressure budget for one Ø 400 mm unit at 0.4 m/s", fontsize=11)
ax4.set_yticks([])
ax4.legend(frameon=True, fontsize=9, loc="upper right",
           bbox_to_anchor=(1.0, 1.35))
ax4.set_xlim(0, left * 1.1)

fig4.savefig("model6_fig4_pressure_budget.png", dpi=200, bbox_inches="tight")
plt.close(fig4)
