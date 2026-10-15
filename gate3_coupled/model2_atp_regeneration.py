"""
Model 2 — ATP Regeneration and Reagent Life
Model 2A: coating liquid alone, showing PPK2 necessity.
Model 2B: reservoir-coupled 7-day run for 10/30/100 L reservoirs.

All concentrations in mM. Time in minutes (2A) or hours/days (2B).
"""

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# --- GS parameters (H195Y, Fisher & Wray 2006) ---

KM_GLU_MM = 30.0
KM_ATP_MM = 7.0
A_GS_UMOL_MIN_MG = 4.7       # specific activity
KM_NH3_MM = 0.4               # Nocek et al. (2008)
KI_GLN_MM = 19.167            # Cheng–Prusoff derived (Model 1)
GLU_ASSAY_MM = 150.0          # used only in Ki derivation, not here

# --- PPK2 parameters (SlPPK2, Sun et al. 2021, Table 2) ---

A_PPK2_UMOL_MIN_MG = 283.0
KM_ADP_MM = 0.61
KM_POLYP_MM = 1.96            # unverified; sensitivity sweep (nishvasa session) shows no effect

# --- Device loading ---

M_KG_M3 = 2.0
ETA_BEAD = 1.00               # Model 3 corrected result at 10 µm
ETA_COAT = 0.985              # Model 4
RHO_RETAINED = 0.5
E_ACTIVE_MG_ML = ETA_BEAD * ETA_COAT * RHO_RETAINED * M_KG_M3  # kg/m^3 ≡ mg/mL

VMAX_GS = A_GS_UMOL_MIN_MG * E_ACTIVE_MG_ML   # mM/min
VMAX_PPK2 = A_PPK2_UMOL_MIN_MG * E_ACTIVE_MG_ML

# --- Reservoir charge ---

GLU_0_MM = 150.0
ATP_0_MM = 3.5
ADP_0_MM = 0.0
POLYP_0_MM = 30.0
GLN_0_MM = 0.0
PI_0_MM = 0.0

# --- Model 2A constants ---

NH3_DISSOLVED_MM = 2.95       # equilibrium with 50 ppm via Henry's law
V_COAT_L = 0.64e-3            # 0.64 L = 0.00064 m^3 -> but we work in mM so volume cancels in rate
PHOSPHATE_INHIB_THRESHOLD_MM = 12.5  # BNU-China 2025


def gs_rate_mM_min(nh3, gln, glu, atp):
    """GS rate in mM/min."""
    km_glu_app = KM_GLU_MM * (1.0 + gln / KI_GLN_MM)
    f_atp = atp / (KM_ATP_MM + atp)
    f_nh3 = nh3 / (KM_NH3_MM + nh3)
    f_glu = glu / (km_glu_app + glu)
    return VMAX_GS * f_atp * f_nh3 * f_glu


def ppk2_rate_mM_min(adp, polyp):
    """PPK2 rate in mM/min."""
    if adp <= 0 or polyp <= 0:
        return 0.0
    f_adp = adp / (KM_ADP_MM + adp)
    f_polyp = polyp / (KM_POLYP_MM + polyp)
    return VMAX_PPK2 * f_adp * f_polyp


# ============================================================
# Model 2A — Coating liquid alone, fixed dissolved NH3
# ============================================================

def ode_2a(t, y, with_ppk2):
    """y = [glu, gln, pi, atp, adp, polyp]. NH3 held at 2.95 mM."""
    glu, gln, pi, atp, adp, polyp = np.maximum(y, 0.0)
    r_gs = gs_rate_mM_min(NH3_DISSOLVED_MM, gln, glu, atp)
    r_ppk2 = ppk2_rate_mM_min(adp, polyp) if with_ppk2 else 0.0
    return [
        -r_gs,                # d[Glu]/dt
        +r_gs,                # d[Gln]/dt
        +r_gs,                # d[Pi]/dt
        -r_gs + r_ppk2,       # d[ATP]/dt
        +r_gs - r_ppk2,       # d[ADP]/dt
        -r_ppk2,              # d[PolyP]/dt
    ]


y0_2a = [GLU_0_MM, GLN_0_MM, PI_0_MM, ATP_0_MM, ADP_0_MM, POLYP_0_MM]
t_2a = (0, 60)  # 60 min
t_eval_2a = np.linspace(0, 60, 2000)

sol_no_ppk2 = solve_ivp(ode_2a, t_2a, y0_2a, args=(False,),
                         method="LSODA", t_eval=t_eval_2a, rtol=1e-8, atol=1e-10)
sol_with_ppk2 = solve_ivp(ode_2a, t_2a, y0_2a, args=(True,),
                           method="LSODA", t_eval=t_eval_2a, rtol=1e-8, atol=1e-10)

# When does ATP drop below 0.5 mM without PPK2?
atp_no = sol_no_ppk2.y[3]
idx_05 = np.argmax(atp_no < 0.5)
if idx_05 > 0:
    print(f"Model 2A: Without PPK2, ATP < 0.5 mM at t = {sol_no_ppk2.t[idx_05]:.1f} min")

# ATP with PPK2
atp_with = sol_with_ppk2.y[3]
print(f"Model 2A: With PPK2, ATP minimum = {atp_with.min():.2f} mM")

# PolyP exhaustion time
polyp_with = sol_with_ppk2.y[5]
idx_polyp = np.argmax(polyp_with < 0.01)
if idx_polyp > 0:
    print(f"Model 2A: PolyP exhausted at t = {sol_with_ppk2.t[idx_polyp]:.0f} min")

# --- Figure 2A: with/without PPK2 comparison ---

fig_2a, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(14, 5), constrained_layout=True)

ax_a.plot(sol_no_ppk2.t, sol_no_ppk2.y[3], color="#D55E00", linewidth=2, label="Without PPK2")
ax_a.plot(sol_with_ppk2.t, sol_with_ppk2.y[3], color="#0072B2", linewidth=2, label="With PPK2")
ax_a.set_xlabel("Time (min)", fontsize=11)
ax_a.set_ylabel("[ATP] (mM)", fontsize=11)
ax_a.set_title("(a) ATP pool", fontsize=12)
ax_a.legend(frameon=True, fontsize=10)
ax_a.set_xlim(0, 60)
ax_a.set_ylim(bottom=0)

ax_b.plot(sol_with_ppk2.t, sol_with_ppk2.y[5], color="#009E73", linewidth=2, label="PolyP")
ax_b.plot(sol_with_ppk2.t, sol_with_ppk2.y[1], color="#CC79A7", linewidth=2, label="Glutamine")
ax_b.plot(sol_with_ppk2.t, sol_with_ppk2.y[2], color="#56B4E9", linewidth=2, label=r"P$_i$")
ax_b.set_xlabel("Time (min)", fontsize=11)
ax_b.set_ylabel("Concentration (mM)", fontsize=11)
ax_b.set_title("(b) Coupled system species", fontsize=12)
ax_b.legend(frameon=True, fontsize=10)
ax_b.set_xlim(0, 60)
ax_b.set_ylim(bottom=0)

fig_2a.savefig("model2a_ppk2_necessity.png", dpi=200, bbox_inches="tight")
plt.close(fig_2a)


# ============================================================
# Model 2B — Reservoir-coupled 7-day simulation
# ============================================================

# Farm parameters: 1,500-bird farm, single unit, batch 80->20 ppm
# Daily NH3 load per unit: 0.307 mol/day
# We treat each day as: shed starts at 80 ppm, unit runs until 20 ppm or 24 h.
# R(C) from model 6: at device loading, approximately 713 µmol/min at 50 ppm.
# For the reservoir-level sim, use the 3-substrate rate law at device E_active.
# The daily NH3 captured = integral of R over the batch time.
# Approximate: each unit captures about 0.307 mol/day (from handoff).

DAILY_LOAD_MOL = 0.307
DAILY_LOAD_MM_PER_L = DAILY_LOAD_MOL * 1000.0  # mmol per day per unit

V_COAT_L_ACTUAL = 0.64  # litres


def simulate_reservoir_7days(v_reservoir_L):
    """Simulate 7 days of batch operation with a given reservoir volume.
    Returns: daily scrub times (h), final composition dict, daily Pi levels.
    """
    v_total_L = v_reservoir_L + V_COAT_L_ACTUAL

    # reservoir composition, shared across days (mM in total volume)
    # initial concentrations diluted into total volume
    glu = GLU_0_MM  # reservoir is charged at these concentrations
    gln = GLN_0_MM
    pi = PI_0_MM
    atp = ATP_0_MM
    polyp = POLYP_0_MM

    scrub_times_h = []
    compositions = []  # (glu, gln, pi, atp, polyp) at end of each day

    for day in range(7):
        # each day: consume DAILY_LOAD_MM_PER_L / v_total_L mmol per litre
        # stoichiometry: 1 mol NH3 -> 1 mol Glu consumed, 1 mol Gln produced,
        #                1 mol Pi produced, 1 mol PolyP bond consumed
        # ATP and ADP cycle (PPK2 keeps ATP buffered, net zero)

        daily_delta_mM = DAILY_LOAD_MM_PER_L / v_total_L  # mM consumed/produced per day

        # check if PolyP can cover this day
        can_scrub = polyp >= daily_delta_mM and glu >= daily_delta_mM

        if not can_scrub:
            scrub_times_h.append(np.nan)
            compositions.append((glu, gln, pi, atp, polyp))
            continue

        glu -= daily_delta_mM
        gln += daily_delta_mM
        pi += daily_delta_mM
        polyp -= daily_delta_mM

        # scrub time depends on current saturation fractions
        km_glu_app = KM_GLU_MM * (1.0 + gln / KI_GLN_MM)
        f_glu = glu / (km_glu_app + glu)
        f_atp = atp / (KM_ATP_MM + atp)
        # relative rate compared to day 0 fresh charge
        km_app_0 = KM_GLU_MM * (1.0 + GLN_0_MM / KI_GLN_MM)
        f_glu_0 = GLU_0_MM / (km_app_0 + GLU_0_MM)
        f_atp_0 = ATP_0_MM / (KM_ATP_MM + ATP_0_MM)
        rate_ratio = (f_glu * f_atp) / (f_glu_0 * f_atp_0)

        # baseline scrub time for this farm at nominal: ~7.5 h (from handoff, 100L day 1)
        baseline_h = 7.5
        scrub_time = baseline_h / rate_ratio

        if scrub_time > 24.0:
            scrub_times_h.append(np.nan)
        else:
            scrub_times_h.append(scrub_time)

        compositions.append((glu, gln, pi, atp, polyp))

    return scrub_times_h, compositions


reservoir_sizes = [10, 30, 100]
reservoir_colors = {"10": "#D55E00", "30": "#E69F00", "100": "#0072B2"}

print("\n--- Model 2B: 7-day reservoir simulation ---")

all_scrub_times = {}
all_compositions = {}
for vr in reservoir_sizes:
    times, comps = simulate_reservoir_7days(vr)
    all_scrub_times[vr] = times
    all_compositions[vr] = comps
    print(f"\nReservoir {vr} L:")
    for d, t in enumerate(times):
        status = f"{t:.1f} h" if not np.isnan(t) else "target not reached"
        print(f"  Day {d+1}: {status}")
    if comps:
        last = comps[-1]
        print(f"  Day 7 composition: Glu={last[0]:.1f}, Gln={last[1]:.1f}, "
              f"Pi={last[2]:.1f}, ATP={last[3]:.1f}, PolyP={last[4]:.1f} mM")


# --- Figure 2B-1: daily scrub time ---

fig_2b1, ax_2b1 = plt.subplots(figsize=(8, 5), constrained_layout=True)
days = np.arange(1, 8)

for vr in reservoir_sizes:
    times = all_scrub_times[vr]
    valid = [(d, t) for d, t in zip(days, times) if not np.isnan(t)]
    failed = [(d, t) for d, t in zip(days, times) if np.isnan(t)]
    c = reservoir_colors[str(vr)]
    if valid:
        vd, vt = zip(*valid)
        ax_2b1.plot(vd, vt, 'o-', color=c, linewidth=2, markersize=6, label=f"{vr} L reservoir")
    if failed:
        for fd in [f[0] for f in failed]:
            ax_2b1.axvline(fd, color=c, linestyle=':', alpha=0.4)
            ax_2b1.plot(fd, 24, 'x', color=c, markersize=10, markeredgewidth=2)

ax_2b1.set_xlabel("Day", fontsize=11)
ax_2b1.set_ylabel("Scrub time (h)", fontsize=11)
ax_2b1.set_title("Daily scrub time — one unit, 1,500-bird farm", fontsize=12)
ax_2b1.legend(frameon=True, fontsize=10, loc="upper left")
ax_2b1.set_xlim(0.5, 7.5)
ax_2b1.set_ylim(0, 26)
ax_2b1.set_xticks(days)
fig_2b1.savefig("model2b_daily_scrub_time.png", dpi=200, bbox_inches="tight")
plt.close(fig_2b1)

# --- Figure 2B-2: reservoir composition for 100 L ---

comps_100 = all_compositions[100]
glu_arr = [c[0] for c in comps_100]
gln_arr = [c[1] for c in comps_100]
pi_arr = [c[2] for c in comps_100]
polyp_arr = [c[4] for c in comps_100]

fig_2b2, ax_2b2 = plt.subplots(figsize=(8, 5), constrained_layout=True)
ax_2b2.plot(days, glu_arr, 'o-', color="#0072B2", linewidth=2, label="Glutamate")
ax_2b2.plot(days, gln_arr, 's-', color="#CC79A7", linewidth=2, label="Glutamine")
ax_2b2.plot(days, pi_arr, '^-', color="#56B4E9", linewidth=2, label=r"P$_i$")
ax_2b2.plot(days, polyp_arr, 'D-', color="#009E73", linewidth=2, label="PolyP")
ax_2b2.axhline(PHOSPHATE_INHIB_THRESHOLD_MM, color="gray", linestyle="--", alpha=0.7,
               label=f"Pi threshold ({PHOSPHATE_INHIB_THRESHOLD_MM} mM)")
ax_2b2.set_xlabel("Day", fontsize=11)
ax_2b2.set_ylabel("Concentration (mM)", fontsize=11)
ax_2b2.set_title("Reservoir composition — 100 L, one unit, 1,500-bird farm", fontsize=12)
ax_2b2.legend(frameon=True, fontsize=10, loc="center right")
ax_2b2.set_xlim(0.5, 7.5)
ax_2b2.set_xticks(days)
fig_2b2.savefig("model2b_reservoir_composition_100L.png", dpi=200, bbox_inches="tight")
plt.close(fig_2b2)

# --- Figure 2B-3: service life vs reservoir volume ---

vr_sweep = np.linspace(5, 200, 200)
v_coat = V_COAT_L_ACTUAL

# PolyP exhaustion: days = V_total * [PolyP]_0 / daily_load_mmol
# where daily_load_mmol = 307 mmol
polyp_life_days = (vr_sweep + v_coat) * POLYP_0_MM / DAILY_LOAD_MM_PER_L

# Phosphate threshold: Pi reaches 12.5 mM after N days
# N * daily_delta_mM = 12.5, daily_delta_mM = 307 / (vr + v_coat)
pi_life_days = PHOSPHATE_INHIB_THRESHOLD_MM * (vr_sweep + v_coat) / DAILY_LOAD_MM_PER_L

ratio = POLYP_0_MM / PHOSPHATE_INHIB_THRESHOLD_MM
print(f"\n[PolyP]_0 / Pi threshold = {ratio:.1f} — PolyP exhaustion is always {ratio:.1f}x "
      f"later than phosphate threshold")

# 7-day line crossings
vr_7day_polyp = 7 * DAILY_LOAD_MM_PER_L / POLYP_0_MM - v_coat
vr_7day_pi = 7 * DAILY_LOAD_MM_PER_L / PHOSPHATE_INHIB_THRESHOLD_MM - v_coat
print(f"For 7-day life: need {vr_7day_polyp:.0f} L (PolyP only) or {vr_7day_pi:.0f} L (Pi constraint)")
print("Phosphate accumulation, not polyphosphate stoichiometry, sets the refill cadence.")

fig_2b3, ax_2b3 = plt.subplots(figsize=(8, 5), constrained_layout=True)
ax_2b3.plot(vr_sweep, polyp_life_days, color="#009E73", linewidth=2, label="PolyP exhaustion")
ax_2b3.plot(vr_sweep, pi_life_days, color="#D55E00", linewidth=2,
            label=f"Pi reaches {PHOSPHATE_INHIB_THRESHOLD_MM} mM")
ax_2b3.axhline(7, color="gray", linestyle="--", alpha=0.6, label="7-day target")
ax_2b3.set_xlabel("Reservoir volume (L)", fontsize=11)
ax_2b3.set_ylabel("Service life (days)", fontsize=11)
ax_2b3.set_title("Service life vs reservoir volume\n"
                  f"(PolyP line is always {ratio:.1f}× later than Pi line)", fontsize=12)
ax_2b3.legend(frameon=True, fontsize=10)
ax_2b3.set_xlim(0, 200)
ax_2b3.set_ylim(bottom=0)
fig_2b3.savefig("model2b_service_life.png", dpi=200, bbox_inches="tight")
plt.close(fig_2b3)

print("\nConclusion: PPK2 is essential — without it GS stops within minutes. "
      "PolyP sets reagent life; phosphate accumulation crosses 12.5 mM 2.4x earlier "
      "than PolyP exhaustion, making it the binding constraint on refill interval.")
