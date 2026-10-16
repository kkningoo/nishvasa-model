#!/usr/bin/env python3
"""
Interface contact stability: identify native inter-chain residue-pair contacts
in a reference frame, then track what fraction of them survive over the
trajectory. Used for Run 3 (Mg-depleted) -- "does the dimer interface stay
together without the catalytic metals?"

Chain identification: do NOT assume a residue-numbering gap marks the chain
boundary -- in this project's tprs, protein residues are numbered
continuously across both chains with NO gap or restart (1..886 straight
through), so the two chains are only distinguishable via MDAnalysis's
`segid`/moltype info exposed by the TPRParser (confirmed present here:
`seg_0_Protein_chain_A` = resid 1-443, `seg_1_Protein_chain_B` = resid
444-886). Check `u.segments` for your own system before assuming a split
point -- a resid-gap heuristic silently picked a nonsense 1-residue "chain A"
the first time this was tried here.

Native contact definition: any inter-chain residue pair with a heavy-atom
minimum distance <= NATIVE_CUTOFF_A in the reference frame (frame 0 of the
corrected trajectory, after running 01_pbc_fix_and_rmsd_rmsf_gyrate.sh).
A pair-residue COM pre-filter (15A) keeps the O(443*443) frame-0 scan fast;
only pairs that pass both the coarse COM filter and the real heavy-atom
distance check become tracked "native pairs".
"""
import MDAnalysis as mda
import MDAnalysis.lib.distances as mdist
import numpy as np
import json
import csv
import warnings
warnings.filterwarnings("ignore")

# ---------------------------------------------------------------------------
# CONFIG
# ---------------------------------------------------------------------------
TOPO = "/home/ayush/md_gs_mgdepleted/md_seg10.tpr"
TRAJ = "/home/ayush/md_gs_mgdepleted/mgdep_fixed_nojump.xtc"
OUTDIR = "/home/ayush/md_gs_mgdepleted/"
OUT_TAG = "Mg_depleted"
CHAINA_SEGID = "seg_0_Protein_chain_A"
CHAINB_SEGID = "seg_1_Protein_chain_B"
NATIVE_CUTOFF_A = 4.5
COM_PREFILTER_A = 15.0
# ---------------------------------------------------------------------------

u = mda.Universe(TOPO, TRAJ)
n_frames = len(u.trajectory)
dt_ns = u.trajectory.dt / 1000.0
print("n_frames:", n_frames, "dt_ns:", dt_ns)

protein = u.select_atoms("protein")
chainA_resids = np.unique(u.select_atoms(f"segid {CHAINA_SEGID}").resids)
chainB_resids = np.unique(u.select_atoms(f"segid {CHAINB_SEGID}").resids)
print("chainA resid range:", chainA_resids.min(), "-", chainA_resids.max(), "n=", len(chainA_resids))
print("chainB resid range:", chainB_resids.min(), "-", chainB_resids.max(), "n=", len(chainB_resids))

u.trajectory[0]
box0 = u.trajectory.ts.dimensions

chainA_atoms = protein.select_atoms("resid " + " ".join(map(str, chainA_resids)))
chainB_atoms = protein.select_atoms("resid " + " ".join(map(str, chainB_resids)))

comA = {r: chainA_atoms.select_atoms(f"resid {r}").center_of_geometry() for r in chainA_resids}
comB = {r: chainB_atoms.select_atoms(f"resid {r}").center_of_geometry() for r in chainB_resids}

native_pairs = []
for ra in chainA_resids:
    for rb in chainB_resids:
        if np.linalg.norm(comA[ra] - comB[rb]) > COM_PREFILTER_A:
            continue
        a1 = chainA_atoms.select_atoms(f"resid {ra}")
        a2 = chainB_atoms.select_atoms(f"resid {rb}")
        if mdist.distance_array(a1.positions, a2.positions, box=box0).min() <= NATIVE_CUTOFF_A:
            native_pairs.append((int(ra), int(rb)))

print("n_native_pairs:", len(native_pairs))

pair_groups = [(chainA_atoms.select_atoms(f"resid {ra}"), chainB_atoms.select_atoms(f"resid {rb}"))
               for ra, rb in native_pairs]

retained_frac = np.zeros(n_frames)
times = np.zeros(n_frames)
per_pair_contact = np.zeros((len(native_pairs), n_frames), dtype=bool)

for i, ts in enumerate(u.trajectory):
    box = ts.dimensions
    times[i] = ts.time / 1000.0
    n_ok = 0
    for j, (a1, a2) in enumerate(pair_groups):
        ok = mdist.distance_array(a1.positions, a2.positions, box=box).min() <= NATIVE_CUTOFF_A
        per_pair_contact[j, i] = ok
        if ok:
            n_ok += 1
    retained_frac[i] = n_ok / len(native_pairs)

csv_path = f"{OUTDIR}interface_native_contact_fraction_{OUT_TAG}.csv"
with open(csv_path, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["time_ns", "fraction_native_contacts_retained"])
    for i in range(n_frames):
        w.writerow([f"{times[i]:.3f}", f"{retained_frac[i]:.4f}"])

first10_mask = times <= 10.0
last10_mask = times >= (times.max() - 10.0)

summary = {
    "n_native_residue_pairs": len(native_pairs),
    "native_cutoff_A": NATIVE_CUTOFF_A,
    "fraction_native_retained_overall_mean": float(retained_frac.mean()),
    "fraction_native_retained_first10ns": float(retained_frac[first10_mask].mean()),
    "fraction_native_retained_last10ns": float(retained_frac[last10_mask].mean()),
    "fraction_native_retained_min": float(retained_frac.min()),
}

protein_resnames = {r: protein.select_atoms(f"resid {r}").resnames[0]
                     for r in list(chainA_resids) + list(chainB_resids)}
pair_details = {}
for j, (ra, rb) in enumerate(native_pairs):
    pair_details[f"A{ra}({protein_resnames[ra]})-B{rb}({protein_resnames[rb]})"] = {
        "fraction_in_contact": float(per_pair_contact[j].mean()),
    }
summary["all_native_pairs"] = pair_details

json_path = f"{OUTDIR}interface_contact_summary_{OUT_TAG}.json"
with open(json_path, "w") as f:
    json.dump(summary, f, indent=2)

print("Wrote", csv_path)
print("Wrote", json_path)
for k in ["fraction_native_retained_overall_mean", "fraction_native_retained_first10ns",
          "fraction_native_retained_last10ns", "fraction_native_retained_min"]:
    print(k, summary[k])
