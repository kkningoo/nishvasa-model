#!/usr/bin/env python3
"""
Contact-lifetime analysis: per-frame distance + binding-event statistics for a
set of named atom-group pairs (e.g. catalytic Mg-substrate coordination,
salt-bridge networks), per chain of a homodimer.

This is a TEMPLATE, not a no-argument tool -- the exact residue/atom selections
below are specific to the GS (glutamine synthetase) dimer construct used in
this project (Run 2: WT vs H195Y inhibited complex, glutamine + phosphate +
2x catalytic Mg per chain). Edit the SELECTIONS block for a different system.

Key gotchas this script already handles, found the hard way during this
project's analysis:

  - Residue numbering is NOT always chain-relative. In this project's tpr
    files, protein residues are numbered CONTINUOUSLY across both chains
    (chain A = resid 1-443, chain B = resid 444-886 -- NOT restarting at 1),
    and chain B's "same local position" as chain A's residue X is NOT
    reliably X+443: verify empirically by scanning a residue range for the
    expected sequence motif (see find_chain_offset() below) rather than
    assuming a fixed offset.

  - The bound ligand (e.g. free glutamine/glutamate substrate) often has the
    SAME residue name as a standard amino acid (GLN, GLU), so a naive
    `resname GLN` selection will also match every real in-chain glutamine --
    filter by resid (or segid/moltype, if your tpr exposes clean segids via
    MDAnalysis -- not guaranteed, GROMACS tprs don't always carry per-chain
    segid info depending on how the topology was built; this script's
    find_* helpers are written to work either way).

  - The catalytic Mg site you want ("n2", closest to the substrate) is not
    necessarily the first MG atom in file order -- there can be 2+ Mg ions
    per chain (n1, n2, differing coordination partners). Identify n2 by
    measuring which Mg is actually close to the ligand atom of interest in a
    single frame, not by assumption (see find_catalytic_mg() below).

  - Run this on a PBC-corrected trajectory (see
    01_pbc_fix_and_rmsd_rmsf_gyrate.sh) -- but note that for SHORT-RANGE,
    INTRA-CHAIN contacts like these, the common chain-vs-chain PBC clustering
    bug usually does NOT corrupt the numbers much (confirmed in this project:
    re-running this analysis before and after the PBC fix gave numbers
    matching to 3 decimal places). It's RMSD/RMSF (whole-complex fits) that
    are PBC-sensitive, not these local distances. Still worth running on the
    corrected trajectory for cleanliness/consistency.
"""
import MDAnalysis as mda
import MDAnalysis.lib.distances as mdist
import numpy as np
import json
import csv
import warnings
warnings.filterwarnings("ignore")

# ---------------------------------------------------------------------------
# CONFIG -- edit per system
# ---------------------------------------------------------------------------
TOPO = "/home/ayush/md_gs_wt_inhib/md_seg20.tpr"       # tpr (for masses/topology)
TRAJ = "/home/ayush/md_gs_wt_inhib/wt_fixed_nojump.xtc"  # PBC-corrected trajectory
OUTDIR = "/home/ayush/md_gs_wt_inhib/"
OUT_TAG = "WT_inhib"  # used in output filenames

# Residue IDs for chain A (adjust per system -- these match the GS dimer's
# own internal numbering, confirmed via diagnostic scans, NOT necessarily the
# paper's residue numbering, which may use a different N-terminal offset).
CHAINA_ARG = 61
CHAINA_GLU = 64
CHAINA_ASP = 52

# cutoffs in Angstrom
MG_LIGAND_CUTOFF = 2.6   # tight: direct coordination bond
SALT_BRIDGE_CUTOFF = 4.0  # looser: electrostatic/H-bond network
# ---------------------------------------------------------------------------


def find_chain_offset(u, probe_range=range(1, 1000), motif=None):
    """Scan CA resnames over probe_range and return the resid where a given
    short sequence motif (e.g. the 6 residues spanning your Arg/Glu/Asp triad)
    re-appears -- i.e. chain B's local numbering offset. Prints candidates;
    does not auto-apply (verify by eye, sequences can repeat)."""
    seq = {}
    for rid in probe_range:
        a = u.select_atoms(f"resid {rid} and name CA")
        if len(a):
            seq[rid] = a.resnames[0]
    if motif:
        motif_len = len(motif)
        rids = sorted(seq)
        for i in range(len(rids) - motif_len):
            window = tuple(seq[rids[i + k]] for k in range(motif_len))
            if window == motif:
                print(f"motif match starting at resid {rids[i]}")
    return seq


def find_catalytic_mg(u, mg_selection, ligand_atom_selection):
    """Given a selection covering multiple Mg ions and ONE ligand atom,
    return the resid of the Mg closest to that atom (the catalytic / n2 Mg)."""
    mg_atoms = u.select_atoms(mg_selection)
    lig = u.select_atoms(ligand_atom_selection)
    assert len(lig) == 1, f"expected exactly 1 ligand atom, got {len(lig)}"
    best = None
    for a in mg_atoms:
        d = np.linalg.norm(a.position - lig.positions[0])
        if best is None or d < best[1]:
            best = (a.resid, d)
    return best  # (resid, distance_at_this_frame)


def mindist(g1, g2, box):
    return mdist.distance_array(g1.positions, g2.positions, box=box).min()


def run_contact_lifetime_analysis(pair_defs, cutoffs, labels):
    """pair_defs: dict name -> (atomgroup1, atomgroup2).
    cutoffs/labels: dict name -> cutoff_A / display name."""
    u = mda.Universe(TOPO, TRAJ)
    n_frames = len(u.trajectory)
    dt_ns = u.trajectory.dt / 1000.0
    print("n_frames:", n_frames, "dt_ns:", dt_ns)

    times = np.zeros(n_frames)
    cols = {k: np.zeros(n_frames) for k in pair_defs}

    for i, ts in enumerate(u.trajectory):
        box = ts.dimensions
        times[i] = ts.time / 1000.0
        for key, (g1, g2) in pair_defs.items():
            if len(g1) == 1 and len(g2) == 1:
                cols[key][i] = mdist.distance_array(g1.positions, g2.positions, box=box)[0, 0]
            else:
                cols[key][i] = mindist(g1, g2, box)

    csv_path = f"{OUTDIR}contact_distances_{OUT_TAG}.csv"
    with open(csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["time_ns"] + list(cols.keys()))
        for i in range(n_frames):
            w.writerow([f"{times[i]:.3f}"] + [f"{cols[k][i]:.3f}" for k in cols])

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
            "name": labels[key],
            "cutoff_A": cutoff,
            "n_frames": n_frames,
            "duration_ns": round(n_frames * dt_ns - dt_ns, 3),
            "fraction_in_contact": float(in_contact.mean()),
            "n_binding_events": len(lifetimes_ns),
            "mean_lifetime_ns": float(np.mean(lifetimes_ns)) if lifetimes_ns else 0.0,
            "median_lifetime_ns": float(np.median(lifetimes_ns)) if lifetimes_ns else 0.0,
            "max_lifetime_ns": float(np.max(lifetimes_ns)) if lifetimes_ns else 0.0,
            "min_dist_A": float(dist.min()),
            "mean_dist_A": float(dist.mean()),
            "max_dist_A": float(dist.max()),
        }

    json_path = f"{OUTDIR}contact_lifetime_summary_{OUT_TAG}.json"
    with open(json_path, "w") as f:
        json.dump(summary, f, indent=2)

    print("Wrote", csv_path)
    print("Wrote", json_path)
    for k, v in summary.items():
        print(k, "frac=%.3f" % v["fraction_in_contact"], "events=", v["n_binding_events"],
              "mean_life_ns=%.3f" % v["mean_lifetime_ns"], "mean_dist_A=%.3f" % v["mean_dist_A"])


if __name__ == "__main__":
    u = mda.Universe(TOPO, TRAJ)

    # Example wiring for the GS dimer (adjust resids/selections for your system):
    mgA_resid, _ = find_catalytic_mg(u, "resname MG", f"resid {CHAINA_ARG} and name CA")  # placeholder probe
    print("NOTE: edit CONFIG and the pair_defs block below for your actual system before relying on this output.")

    # pair_defs = {
    #     "chainA_Mg_n2_Gln": (u.select_atoms(f"resid {mgA_resid} and resname MG"),
    #                          u.select_atoms("resid <ligand_resid> and name OE1")),
    #     "chainA_Arg_Glu": (u.select_atoms(f"resid {CHAINA_ARG} and name NE NH1 NH2"),
    #                        u.select_atoms(f"resid {CHAINA_GLU} and name OE1 OE2")),
    #     "chainA_Arg_Asp": (u.select_atoms(f"resid {CHAINA_ARG} and name NE NH1 NH2"),
    #                        u.select_atoms(f"resid {CHAINA_ASP} and name OD1 OD2")),
    #     # ... mirror for chain B with its own (empirically-found) resids ...
    # }
    # cutoffs = {"chainA_Mg_n2_Gln": MG_LIGAND_CUTOFF, "chainA_Arg_Glu": SALT_BRIDGE_CUTOFF, "chainA_Arg_Asp": SALT_BRIDGE_CUTOFF}
    # labels = {"chainA_Mg_n2_Gln": "chainA_Mg(n2)-Gln(OE1)", "chainA_Arg_Glu": "chainA_Arg-Glu", "chainA_Arg_Asp": "chainA_Arg-Asp"}
    # run_contact_lifetime_analysis(pair_defs, cutoffs, labels)
