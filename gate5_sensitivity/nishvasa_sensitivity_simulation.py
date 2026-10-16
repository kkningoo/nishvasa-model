#!/usr/bin/env python3
"""
NISHVASA cartridge model: sensitivity and uncertainty simulation (single file)
==============================================================================

Contents
  1. Model engine (docx Models 3-7): kinetics, Thiele moduli, gas film, plug flow, batch time
  2. Verification against the modelling document and independent solvers
  3. Part A: local elasticity/tornado, Morris screening, Sobol indices (+ convergence replicate)
  4. Part B: Monte Carlo (Latin hypercube) propagation, success probabilities, units for 90 % confidence

Usage
  python nishvasa_sensitivity_simulation.py verify     # baseline check vs document (+ LSODA cross-check)
  python nishvasa_sensitivity_simulation.py partA      # writes sens1_results.pkl
  python nishvasa_sensitivity_simulation.py partB      # writes mc_results.pkl and S2_deployment_table.csv
  python nishvasa_sensitivity_simulation.py all        # verify + partA + partB

Requirements: numpy, scipy, SALib (pip install SALib)

Units: SI internally (mol, m, s). Gas NH3 in ppm, 1 ppm = 4.09e-5 mol/m3 (25 C).
"""
import sys
import time
import csv
import pickle
import numpy as np
from scipy.stats import qmc, norm, triang, truncnorm

# =====================================================================================
# 1. MODEL ENGINE
# =====================================================================================
PPM = 4.09e-5          # mol/m3 per ppm
Z_STACK = 0.280        # m, stack depth (14 trays x 20 mm)
PACK = 0.92            # hexagon-in-circle packing fraction
D_AIR = 2.3e-5         # m2/s, NH3 diffusivity in air
L_BEAD = 10e-6         # m, selected bead radius (Model 3)
D_ATP_FREE = 3.7e-10   # m2/s, free ATP diffusivity in water

# Baseline parameter set (modelling document, Models 1-6)
BASE = dict(
    rho=0.5,                         # activity retention
    ATP=3.5, Glu=150.0,              # mM, reservoir levels
    M=2.0,                           # kg/m3 enzyme loading
    a=4.7,                           # umol/min/mg specific activity (GS H195Y)
    KmATP=7.0, KmNH3=0.4, KmGlu=30.0,  # mM
    Ki=115.0 / (1 + 150.0 / 30.0),   # mM, Cheng-Prusoff from IC50 = 115 mM at 150 mM Glu
    Gln=0.0,                         # mM, glutamine carried in reservoir
    Lcoat=24e-6, r=1.5e-3, tw=0.8e-3,  # m
    Dgel=0.5 * D_ATP_FREE,           # m2/s
    Hcc=1443.0, u=0.4, Sh=3.0,
    C0=80.0,                         # ppm, shed start level
    D=0.400,                         # m, disc diameter
)

_GL_X, _GL_W = np.polynomial.legendre.leggauss(8)   # 8-point Gauss-Legendre nodes/weights


def eta_sphere(phi):
    """Effectiveness factor of a sphere, (3/phi)(coth phi - 1/phi), with a series for small phi."""
    phi = np.asarray(phi, dtype=float)
    out = np.empty_like(phi)
    small = phi < 0.05
    ps = phi[small]
    out[small] = 1.0 - ps ** 2 / 15.0 + 2.0 * ps ** 4 / 315.0
    pl = phi[~small]
    out[~small] = (3.0 / pl) * (1.0 / np.tanh(pl) - 1.0 / pl)
    return out


def derived(p, eta_coat_fixed=None, eta_bead_fixed=None):
    """Derived kinetic and geometric quantities. p: dict of scalars or equal-length arrays
    (missing keys fall back to BASE)."""
    g = lambda k: np.asarray(p.get(k, BASE[k]), dtype=float)
    rho, ATP, Glu, M, a = g("rho"), g("ATP"), g("Glu"), g("M"), g("a")
    KmATP, KmNH3, Ki, Gln = g("KmATP"), g("KmNH3"), g("Ki"), g("Gln")
    Lc, r, Dg, Hcc, u, Sh, tw, KmGlu, D = (g("Lcoat"), g("r"), g("Dgel"), g("Hcc"),
                                            g("u"), g("Sh"), g("tw"), g("KmGlu"), g("D"))
    a_si = a / 60.0                                    # umol/min/mg -> mol/s/kg
    aM = a_si * M                                      # mol s-1 m-3, intrinsic volumetric rate
    k2 = aM / ((KmATP + ATP) * Dg)                     # Thiele modulus squared per length squared (ATP basis)
    phi_c = Lc * np.sqrt(k2)
    phi_b = L_BEAD * np.sqrt(k2)
    eta_c = (np.tanh(phi_c) / phi_c) if eta_coat_fixed is None else np.full_like(phi_c, eta_coat_fixed)
    eta_b = (eta_sphere(np.atleast_1d(phi_b)).reshape(np.shape(phi_b)) if eta_bead_fixed is None
             else np.full_like(phi_b, eta_bead_fixed))
    fGlu = Glu / (KmGlu * (1.0 + Gln / Ki) + Glu)      # competitive inhibition by glutamine
    fATP = ATP / (KmATP + ATP)
    Vmax = eta_b * eta_c * rho * aM * fGlu * fATP      # mol s-1 m-3, effective maximum rate
    # geometry
    ro, rout = r - Lc, r + tw / 2.0
    Aopen = 2 * np.sqrt(3) * ro ** 2
    Acell = 2 * np.sqrt(3) * rout ** 2
    P = 4 * np.sqrt(3) * ro
    Dh = 2 * ro
    eps = Aopen / Acell
    face = np.pi * D ** 2 / 4.0
    Ncell = PACK * face / Acell
    Q = u * face
    uch = u / (PACK * eps)
    kg = Sh * D_AIR / Dh
    return dict(Vmax=Vmax, eta_c=eta_c, eta_b=eta_b, fGlu=fGlu, fATP=fATP, Aopen=Aopen, P=P, Dh=Dh,
                eps=eps, Ncell=Ncell, Q=Q, uch=uch, kg=kg, Lc=Lc, Hcc=Hcc, KmNH3=KmNH3, face=face)


def surface_conc(C, d):
    """Liquid-phase NH3 at the wall, c = Hcc*Cs (mol/m3), from the film balance
    kg (C - c/H) = V' c/(Km + c). Closed-form positive root of the quadratic. C in mol/m3."""
    Vp = d["Vmax"] * d["Lc"]
    H, Km, kg = d["Hcc"], d["KmNH3"], d["kg"]
    b = H * C - Km - Vp * H / kg
    return 0.5 * (b + np.sqrt(b * b + 4.0 * H * C * Km))


def wall_flux(C, d):
    c = surface_conc(C, d)
    return d["Vmax"] * d["Lc"] * c / (d["KmNH3"] + c)


def outlet(Cin, d, nstep=24):
    """Plug flow along one channel, fixed-step RK4. Cin and result in mol/m3."""
    h = Z_STACK / nstep
    kappa = d["P"] / (d["uch"] * d["Aopen"])
    f = lambda C: -kappa * wall_flux(C, d)
    C = np.array(Cin, dtype=float, copy=True)
    for _ in range(nstep):
        k1 = f(C)
        k2 = f(C + 0.5 * h * k1)
        k3 = f(C + 0.5 * h * k2)
        k4 = f(C + h * k3)
        C = C + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6.0
    return C


def removal_rate(C_ppm, d, nstep=24):
    """R = Q (Cin - Cout) for one unit, mol/s."""
    Cin = np.asarray(C_ppm, dtype=float) * PPM
    return d["Q"] * (Cin - outlet(Cin, d, nstep))


def kinetic_limit(C_ppm, d):
    Cin = np.asarray(C_ppm, dtype=float) * PPM
    return wall_flux(Cin, d) * d["Ncell"] * d["P"] * Z_STACK


def inv_rate_integral(C0_ppm, d, Cf_ppm=20.0, nstep=24):
    """I = integral from Cf to C0 of dC / R(C), by 8-point Gauss-Legendre.
    Batch time with n identical units: t = Vshed * I / n."""
    C0 = np.asarray(C0_ppm, dtype=float) * PPM
    Cf = Cf_ppm * PPM
    half = 0.5 * (C0 - Cf)
    mid = 0.5 * (C0 + Cf)
    tot = 0.0
    for x, w in zip(_GL_X, _GL_W):
        Cn = mid + half * x
        R = d["Q"] * (Cn - outlet(Cn, d, nstep))
        tot = tot + w / R
    return half * tot


def batch_time_h(C0_ppm, Vshed, n_units, d, Cf_ppm=20.0, nstep=24):
    """Hours to bring a well-mixed shed (Vshed m3) from C0 to Cf ppm with n recirculating units."""
    return Vshed * inv_rate_integral(C0_ppm, d, Cf_ppm, nstep) / n_units / 3600.0


# =====================================================================================
# 2. VERIFICATION
# =====================================================================================
def verify():
    from scipy.integrate import solve_ivp
    d = derived(dict(BASE))
    d98 = derived(dict(BASE), eta_coat_fixed=0.98, eta_bead_fixed=1.0)   # document rounding
    print("Vmax_eff (computed eta) = %.5f, eta_coat = %.4f, eta_bead = %.5f" %
          (float(d["Vmax"]), float(d["eta_c"]), float(d["eta_b"])))
    print("Vmax_eff (document eta 0.98) = %.5f   [document: 0.021]" % float(d98["Vmax"]))
    for C in (20, 50, 100):
        print("R(%3d ppm) = %.1f umol/min   [document: 585 / 715 / 765]" %
              (C, float(removal_rate(C, d)) * 6e7))
    farms = [("500 birds, 1 unit", 150, 60, 1), ("1500 birds, 3 units", 375, 80, 3), ("3000 birds, 6 units", 750, 100, 6)]
    for lab, V, C0, n in farms:
        print("batch time %-22s = %.2f h   [document: 6.0 / 7.3 / 9.6]" % (lab, float(batch_time_h(C0, V, n, d))))
    # independent LSODA cross-check of the channel integration
    kappa = float(d["P"] / (d["uch"] * d["Aopen"]))
    for C in (20, 50, 100):
        s = solve_ivp(lambda x, c: -kappa * wall_flux(c, d), [0, Z_STACK], [C * PPM],
                      rtol=1e-11, atol=1e-16, method="LSODA")
        R_ref = float(d["Q"]) * (C * PPM - s.y[0, -1])
        print("RK4 vs LSODA at %3d ppm: relative error %.2e" % (C, abs(float(removal_rate(C, d)) / R_ref - 1)))


# =====================================================================================
# 3. PART A: LOCAL, MORRIS, SOBOL
# =====================================================================================
V_SHED_A, N_UNITS_A = 375.0, 3        # 1,500-bird farm, 3 units
CF_PPM = 20.0

# name, low, high, scale ("lin" or "log"). Sampled on the unit cube and mapped to these ranges.
PAR = [
    ("rho",   0.3,   0.7,    "lin"),
    ("ATP",   1.75,  7.0,    "lin"),
    ("Glu",   100.0, 200.0,  "lin"),
    ("M",     1.0,   4.0,    "lin"),
    ("a",     4.7 * 0.8, 4.7 * 1.2, "lin"),
    ("KmATP", 5.0,   9.0,    "lin"),
    ("KmNH3", 0.2,   0.8,    "lin"),
    ("Ki",    BASE["Ki"] * 0.7, BASE["Ki"] * 1.3, "lin"),
    ("Gln",   0.0,   21.0,   "lin"),
    ("Lcoat", 20e-6, 50e-6,  "lin"),          # lower limit: 10 um bead must fit (L_bead <= L_coat/2)
    ("r",     1.3e-3, 1.7e-3, "lin"),
    ("Dgel",  0.3 * D_ATP_FREE, 0.8 * D_ATP_FREE, "log"),
    ("Hcc",   1443 * 0.9, 1443 * 1.1, "lin"),
    ("u",     0.2,   0.6,    "lin"),
    ("Sh",    2.5,   4.0,    "lin"),
    ("C0",    60.0,  100.0,  "lin"),
]
NAMES = [p[0] for p in PAR]


def to_phys(U, names):
    """Map unit-cube samples U (n x len(names)) to physical parameter arrays."""
    out = {}
    for j, n in enumerate(names):
        _, lo, hi, sc = PAR[NAMES.index(n)]
        out[n] = (lo + U[:, j] * (hi - lo) if sc == "lin"
                  else np.exp(np.log(lo) + U[:, j] * (np.log(hi) - np.log(lo))))
    return out


def model_AB(p, chunk=4000):
    """Evaluate both outputs for a dict of parameter arrays.
    Returns R at 50 ppm (umol/min, one unit) and batch time (h) for the 1,500-bird farm with 3 units."""
    n = len(next(iter(p.values())))
    R50 = np.empty(n)
    T = np.empty(n)
    for s in range(0, n, chunk):
        sl = slice(s, s + chunk)
        q = {k: v[sl] for k, v in p.items()}
        d = derived(q)
        R50[sl] = removal_rate(50.0, d) * 6e7
        C0 = q.get("C0", np.full(len(next(iter(q.values()))), BASE["C0"]))
        T[sl] = batch_time_h(C0, V_SHED_A, N_UNITS_A, d, CF_PPM)
    return R50, T


def run_part_A(out="sens1_results.pkl"):
    from SALib.sample import sobol as sobol_sample
    from SALib.sample import morris as morris_sample
    from SALib.analyze import sobol as sobol_an
    from SALib.analyze import morris as morris_an

    res, t0 = {}, time.time()

    # ---------- local tornado (+/-20 %) and elasticities (central difference in logs, +/-1 %)
    base = {k: np.array([BASE[k]]) for k in NAMES}
    R0, T0 = model_AB(base)
    res["base"] = dict(R=float(R0[0]), T=float(T0[0]))
    tor = {}
    for n in NAMES:
        if n in ("Gln", "Ki"):          # inactive at the fresh-reservoir baseline (Gln = 0)
            continue
        row = {}
        for tag, f in (("lo", 0.8), ("hi", 1.2), ("elo", 0.99), ("ehi", 1.01)):
            q = {k: v.copy() for k, v in base.items()}
            q[n] = q[n] * f
            R, T = model_AB(q)
            row[tag] = (float(R[0]), float(T[0]))
        dl = np.log(1.01) - np.log(0.99)
        row["elast"] = (float((np.log(row["ehi"][0]) - np.log(row["elo"][0])) / dl),
                        float((np.log(row["ehi"][1]) - np.log(row["elo"][1])) / dl))
        tor[n] = row
    res["tornado"] = tor
    print("tornado done (%.1f s); baseline R50 = %.1f umol/min, T = %.2f h" % (time.time() - t0, R0[0], T0[0]))

    # ---------- Morris screening: r = 50 trajectories, p = 4 levels, step 2/3 on the unit cube
    def run_morris(names, output, r=50, seed=11):
        prob = dict(num_vars=len(names), names=names, bounds=[[0, 1]] * len(names))
        X = morris_sample.sample(prob, r, num_levels=4, seed=seed)
        R, T = model_AB(to_phys(X, names))
        Y = R if output == "R" else np.log(T)
        a = morris_an.analyze(prob, X, Y, conf_level=0.95, num_levels=4, seed=seed)
        return {k: np.asarray(a[k]) for k in ("mu", "mu_star", "sigma", "mu_star_conf")}, names

    names_R = [n for n in NAMES if n != "C0"]          # C0 is not an input of R at fixed 50 ppm
    res["morris_R"], res["morris_R_names"] = run_morris(names_R, "R")
    res["morris_T"], res["morris_T_names"] = run_morris(NAMES, "T")
    print("morris done (%.1f s)" % (time.time() - t0))

    # ---------- Sobol: Saltelli sampling (N(k+2) runs), Saltelli-2010 S1 and Jansen ST, 1000 bootstraps
    def run_sobol(names, N, seed):
        prob = dict(num_vars=len(names), names=names, bounds=[[0, 1]] * len(names))
        X = sobol_sample.sample(prob, N, calc_second_order=False, scramble=True, seed=seed)
        R, T = model_AB(to_phys(X, names))
        out_ = {}
        for lab, Y in (("R", R), ("lnT", np.log(T)), ("T", T)):
            a = sobol_an.analyze(prob, Y, calc_second_order=False, num_resamples=1000,
                                 conf_level=0.95, seed=seed)
            out_[lab] = {k: np.asarray(a[k]) for k in ("S1", "S1_conf", "ST", "ST_conf")}
            out_[lab]["Ymean"], out_[lab]["Yvar"] = float(Y.mean()), float(Y.var())
        out_["n_eval"] = len(X)
        return out_

    def sobol_pair(N, seed):
        oR = run_sobol(names_R, N, seed)               # output R (no C0)
        oT = run_sobol(NAMES, N, seed)                 # output ln t (all 16 inputs)
        return {"R": oR["R"], "n_R": oR["n_eval"], "lnT": oT["lnT"], "T": oT["T"], "n_T": oT["n_eval"]}

    res["sobol_names_R"] = names_R
    res["sobol_N2048"] = sobol_pair(2048, seed=7)
    print("sobol N=2048 done (%.1f s)" % (time.time() - t0))
    res["sobol_N1024"] = sobol_pair(1024, seed=101)    # independent scramble = convergence check
    print("sobol N=1024 done (%.1f s)" % (time.time() - t0))

    # ---------- scatter data: independent Latin hypercube over the same ranges
    U = qmc.LatinHypercube(d=len(NAMES), seed=3).random(4000)
    R, T = model_AB(to_phys(U, NAMES))
    res["scatter"] = dict(X=to_phys(U, NAMES), R=R, T=T)
    res["time"] = time.time() - t0
    with open(out, "wb") as fh:
        pickle.dump(res, fh)

    # short console summary
    s = res["sobol_N2048"]["lnT"]
    print("\nSobol (ln t), N = 2048:  sum S1 = %.3f, sum ST = %.3f" % (s["S1"].sum(), s["ST"].sum()))
    for i in np.argsort(-s["ST"])[:8]:
        print("  %-6s S1 = %.3f  ST = %.3f +/- %.3f" % (NAMES[i], s["S1"][i], s["ST"][i], s["ST_conf"][i]))
    print("saved", out)
    return res


# =====================================================================================
# 4. PART B: MONTE CARLO UNCERTAINTY PROPAGATION
# =====================================================================================
N_MC = 10000
FARMS = {500: dict(V=150.0, C0=60.0, n_rec=1),
         1500: dict(V=375.0, C0=80.0, n_rec=3),
         3000: dict(V=750.0, C0=100.0, n_rec=6)}
C0_HALF_RANGE = 15.0        # ppm, uniform around the document start level
TLIM = 10.0                 # h, acceptance target
SCEN = {"fresh": dict(Glu=150.0, Gln=0.0),     # day 1 of a new reservoir
        "week":  dict(Glu=129.0, Gln=21.0)}    # day 7 of a 100 L reservoir (Model 2B)
KI0 = BASE["Ki"]
MC_COLS = ["rho", "ATP", "a", "KmATP", "Ki", "KmNH3", "Dgel", "Lcoat", "Hcc", "Sh", "C0"]


def mc_draw(seed=2026, n=N_MC):
    """Latin hypercube draw mapped through inverse CDFs. Returns parameter dict and U(0,1) for C0."""
    U = qmc.LatinHypercube(d=len(MC_COLS), seed=seed).random(n)
    u = {c: U[:, i] for i, c in enumerate(MC_COLS)}
    p = dict(
        rho=triang(c=0.5, loc=0.3, scale=0.4).ppf(u["rho"]),                   # triangular(0.3, 0.5, 0.7)
        ATP=2.8 + (4.2 - 2.8) * u["ATP"],                                        # uniform, +/-20 % of 3.5 mM
        a=4.7 * np.exp(0.15 * norm.ppf(u["a"])),                                 # log-normal, sigma_ln 0.15
        KmATP=7.0 * np.exp(0.15 * norm.ppf(u["KmATP"])),                         # log-normal, sigma_ln 0.15
        Ki=KI0 * np.exp(0.25 * norm.ppf(u["Ki"])),                               # log-normal, sigma_ln 0.25
        KmNH3=0.4 * np.exp(0.30 * norm.ppf(u["KmNH3"])),                         # log-normal, sigma_ln 0.30
        Dgel=np.exp(np.log(0.3 * D_ATP_FREE) +
                    u["Dgel"] * (np.log(0.8 * D_ATP_FREE) - np.log(0.3 * D_ATP_FREE))),   # log-uniform
        Lcoat=truncnorm(-1.6, 2.4, loc=24.0, scale=2.5).ppf(u["Lcoat"]) * 1e-6,  # N(24, 2.5) um on [20, 30]
        Hcc=1443.0 * (0.9 + 0.2 * u["Hcc"]),
        Sh=2.5 + 1.5 * u["Sh"],
    )
    return p, u["C0"]


def mc_run(seed=2026):
    p, uC0 = mc_draw(seed)
    out = {"params": p, "uC0": uC0}
    for sname, sc in SCEN.items():
        q = dict(p)
        q["Glu"] = np.full(N_MC, sc["Glu"])
        q["Gln"] = np.full(N_MC, sc["Gln"])
        d = derived(q)
        for farm, f in FARMS.items():
            C0 = f["C0"] - C0_HALF_RANGE + 2 * C0_HALF_RANGE * uC0
            T1 = batch_time_h(C0, f["V"], 1, d, 20.0)            # hours with ONE unit; t_n = T1 / n exactly
            T1n = batch_time_h(np.full(N_MC, f["C0"]), f["V"], 1, d, 20.0)   # start level fixed at nominal
            out[(sname, farm)] = dict(T1=T1, C0=C0, T1_nomC0=T1n)
        out[(sname, "R50")] = removal_rate(50.0, d) * 6e7        # umol/min, one unit at 50 ppm
    # nominal (point-estimate) reference
    nom = {}
    for sname, sc in SCEN.items():
        q = {k: np.array([BASE[k]]) for k in ("rho", "ATP", "a", "KmATP", "Ki", "KmNH3", "Dgel", "Lcoat", "Hcc", "Sh")}
        q["Glu"], q["Gln"] = np.array([sc["Glu"]]), np.array([sc["Gln"]])
        d = derived(q)
        for farm, f in FARMS.items():
            nom[(sname, farm)] = float(batch_time_h(np.array([f["C0"]]), f["V"], 1, d, 20.0)[0])
    out["nominal_T1"] = nom
    return out


def wilson(k, n, z=1.96):
    """95 % Wilson score interval for a binomial proportion."""
    ph = k / n
    den = 1 + z * z / n
    c = (ph + z * z / (2 * n)) / den
    h = z * np.sqrt(ph * (1 - ph) / n + z * z / (4 * n * n)) / den
    return c - h, c + h


def run_part_B(out="mc_results.pkl", csv_out="S2_deployment_table.csv", nmax=12):
    t0 = time.time()
    res = mc_run()
    summ = {}
    for sname in SCEN:
        for farm, f in FARMS.items():
            T1 = res[(sname, farm)]["T1"]
            rows = []
            for n in range(1, nmax + 1):
                t = T1 / n
                k = int((t <= TLIM).sum())
                lo, hi = wilson(k, N_MC)
                rows.append(dict(n=n, P=k / N_MC, P_lo=lo, P_hi=hi,
                                 P_nomC0=float(((res[(sname, farm)]["T1_nomC0"] / n) <= TLIM).mean()),
                                 med=float(np.median(t)), p5=float(np.quantile(t, 0.05)),
                                 p95=float(np.quantile(t, 0.95)), nominal=res["nominal_T1"][(sname, farm)] / n))
            n90 = next((r["n"] for r in rows if r["P"] >= 0.9), None)          # smallest n with P >= 0.9
            n_i = np.ceil(T1 / TLIM).astype(int)                               # units each sample would need
            n90_q = int(np.ceil(np.quantile(T1 / TLIM, 0.9)))                  # equivalent quantile form
            summ[(sname, farm)] = dict(rows=rows, n90=n90, n_i_q90=n90_q, n_i=n_i)
            r = rows[f["n_rec"] - 1]
            print("%-5s %4d birds: P(design %d units) = %.3f | median %.2f h (%.2f to %.2f) | n90 = %s (quantile form %d)" %
                  (sname, farm, f["n_rec"], r["P"], r["med"], r["p5"], r["p95"], n90, n90_q))
    res["summary"] = summ
    with open(out, "wb") as fh:
        pickle.dump(res, fh)

    with open(csv_out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["farm_birds", "reservoir", "units_design", "P_success_design_%", "P_success_design_fixedC0_%",
                    "median_h", "p5_h", "p95_h", "nominal_h", "units_for_90pct", "P_at_n90_%", "n90_CI95_lo_%", "n90_CI95_hi_%"])
        for sc in ("fresh", "week"):
            for farm, f in FARMS.items():
                rows = summ[(sc, farm)]["rows"]
                r, n90 = rows[f["n_rec"] - 1], summ[(sc, farm)]["n90"]
                r9 = rows[n90 - 1] if n90 else None
                w.writerow([farm, sc, f["n_rec"], f"{100 * r['P']:.1f}", f"{100 * r['P_nomC0']:.1f}",
                            f"{r['med']:.2f}", f"{r['p5']:.2f}", f"{r['p95']:.2f}", f"{r['nominal']:.2f}",
                            n90, f"{100 * r9['P']:.1f}" if r9 else "", f"{100 * r9['P_lo']:.1f}" if r9 else "",
                            f"{100 * r9['P_hi']:.1f}" if r9 else ""])
    print("saved %s and %s (%.1f s)" % (out, csv_out, time.time() - t0))
    return res


# =====================================================================================
if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "all"
    if mode in ("verify", "all"):
        verify()
    if mode in ("partA", "all"):
        run_part_A()
    if mode in ("partB", "all"):
        run_part_B()
    if mode not in ("verify", "partA", "partB", "all"):
        print(__doc__)
