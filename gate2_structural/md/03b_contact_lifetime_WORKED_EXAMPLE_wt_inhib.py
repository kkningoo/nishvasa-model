#!/usr/bin/env python3
"""
Concrete, actually-ran version of 03_contact_lifetime_analysis.py for Run 2
(WT inhibited, glutamine + phosphate + 2 catalytic Mg per chain). Kept as a
worked example since it's easier to adapt a real filled-in script than a
template with placeholders.

Residue IDs below were found empirically for this system (continuous
numbering, chain A = resid 1-443, chain B = resid 444-886, offset +443 not
+444 -- confirmed by scanning for the Asp-Gly-Ser-Ser-Ile-Glu...Arg...Glu-Glu
sequence motif around resid 489-510). Catalytic Mg(n2) resids (890 for chain
A, 892 for chain B) were found by measuring which of the 4 Mg ions in the
system sits ~2A from each chain's bound glutamine OE1 (the other two, 891/893,
are the n1 site, further away, near the phosphate).
"""
import MDAnalysis as mda
import numpy as np
import json
import csv
import warnings
warnings.filterwarnings("ignore")

topo = "/home/ayush/md_gs_wt_inhib/md_seg20.tpr"
traj = "/home/ayush/md_gs_wt_inhib/wt_fixed_nojump.xtc"
outdir = "/home/ayush/md_gs_wt_inhib/"

u = mda.Universe(topo, traj)
n_frames = len(u.trajectory)
dt_ns = u.trajectory.dt / 1000.0
print("n_frames:", n_frames, "dt_ns:", dt_ns)

mgA = u.select_atoms("resid 890 and resname MG")
glnA = u.select_atoms("resid 887 and name OE1")
arg61A = u.select_atoms("resid 61 and name NE NH1 NH2")
glu64A = u.select_atoms("resid 64 and name OE1 OE2")
asp52A = u.select_atoms("resid 52 and name OD1 OD2")

mgB = u.select_atoms("resid 892 and resname MG")
glnB = u.select_atoms("resid 888 and name OE1")
arg61B = u.select_atoms("resid 504 and name NE NH1 NH2")
glu64B = u.select_atoms("resid 507 and name OE1 OE2")
asp52B = u.select_atoms("resid 495 and name OD1 OD2")

for name, g, exp in [("mgA", mgA, 1), ("glnA", glnA, 1), ("arg61A", arg61A, 3),
                      ("glu64A", glu64A, 2), ("asp52A", asp52A, 2),
                      ("mgB", mgB, 1), ("glnB", glnB, 1), ("arg61B", arg61B, 3),
                      ("glu64B", glu64B, 2), ("asp52B", asp52B, 2)]:
    assert len(g) == exp, f"{name} expected {exp} got {len(g)}"


def mindist(g1, g2, box):
    d = mda.lib.distances.distance_array(g1.positions, g2.positions, box=box)
    return d.min()


times = np.zeros(n_frames)
cols = {k: np.zeros(n_frames) for k in
        ["chainA_Mg_n2_Gln", "chainA_Arg61_Glu64", "chainA_Arg61_Asp52",
         "chainB_Mg_n2_Gln", "chainB_Arg61_Glu64", "chainB_Arg61_Asp52"]}

for i, ts in enumerate(u.trajectory):
    box = ts.dimensions
    times[i] = ts.time / 1000.0
    cols["chainA_Mg_n2_Gln"][i] = mda.lib.distances.distance_array(mgA.positions, glnA.positions, box=box)[0, 0]
    cols["chainA_Arg61_Glu64"][i] = mindist(arg61A, glu64A, box)
    cols["chainA_Arg61_Asp52"][i] = mindist(arg61A, asp52A, box)
    cols["chainB_Mg_n2_Gln"][i] = mda.lib.distances.distance_array(mgB.positions, glnB.positions, box=box)[0, 0]
    cols["chainB_Arg61_Glu64"][i] = mindist(arg61B, glu64B, box)
    cols["chainB_Arg61_Asp52"][i] = mindist(arg61B, asp52B, box)

csv_path = outdir + "contact_distances_WT_inhib.csv"
with open(csv_path, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["time_ns"] + list(cols.keys()))
    for i in range(n_frames):
        w.writerow([f"{times[i]:.3f}"] + [f"{cols[k][i]:.3f}" for k in cols])

cutoffs = {"chainA_Mg_n2_Gln": 2.6, "chainB_Mg_n2_Gln": 2.6,
           "chainA_Arg61_Glu64": 4.0, "chainB_Arg61_Glu64": 4.0,
           "chainA_Arg61_Asp52": 4.0, "chainB_Arg61_Asp52": 4.0}
labels = {"chainA_Mg_n2_Gln": "chainA_Mg(n2)-Gln(OE1)", "chainB_Mg_n2_Gln": "chainB_Mg(n2)-Gln(OE1)",
          "chainA_Arg61_Glu64": "chainA_Arg61-Glu64", "chainB_Arg61_Glu64": "chainB_Arg61-Glu64",
          "chainA_Arg61_Asp52": "chainA_Arg61-Asp52", "chainB_Arg61_Asp52": "chainB_Arg61-Asp52"}

summary = {}
for key, dist in cols.items():
    cutoff = cutoffs[key]
    in_contact = dist <= cutoff
    events, run_len = [], 0
    for v in in_contact:
        if v:
            run_len += 1
        else:
            if run_len > 0:
                events.append(run_len)
            run_len = 0
    if run_len > 0:
        events.append(run_len)
    lifetimes_ns = [e * dt_ns for e in events]
    summary[key] = {
        "name": labels[key], "cutoff_A": cutoff, "n_frames": n_frames,
        "duration_ns": round(n_frames * dt_ns - dt_ns, 3),
        "fraction_in_contact": float(in_contact.mean()),
        "n_binding_events": len(lifetimes_ns),
        "mean_lifetime_ns": float(np.mean(lifetimes_ns)) if lifetimes_ns else 0.0,
        "median_lifetime_ns": float(np.median(lifetimes_ns)) if lifetimes_ns else 0.0,
        "max_lifetime_ns": float(np.max(lifetimes_ns)) if lifetimes_ns else 0.0,
        "min_dist_A": float(dist.min()), "mean_dist_A": float(dist.mean()), "max_dist_A": float(dist.max()),
    }

json_path = outdir + "contact_lifetime_summary_WT_inhib.json"
with open(json_path, "w") as f:
    json.dump(summary, f, indent=2)

print("Wrote", csv_path)
print("Wrote", json_path)
for k, v in summary.items():
    print(k, "frac=%.3f" % v["fraction_in_contact"], "events=", v["n_binding_events"],
          "mean_life_ns=%.3f" % v["mean_lifetime_ns"], "mean_dist_A=%.3f" % v["mean_dist_A"])
