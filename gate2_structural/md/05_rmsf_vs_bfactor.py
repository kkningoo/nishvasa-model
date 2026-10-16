#!/usr/bin/env python3
"""
Compare simulated per-residue RMSF against crystallographic B-factors
(validation run, Run 1: WT holo vs 4LNK).

Conversion: B_sim = (8*pi^2/3) * RMSF_Angstrom^2 (the standard
isotropic-displacement relation; RMSF must be converted nm -> Angstrom first).

Gotchas found while building this for the GS dimer:
  - gmx rmsf's printed residue-number column is NOT globally unique -- it
    restarts at 1 for chain B (same underlying quirk documented in
    01_pbc_fix_and_rmsd_rmsf_gyrate.sh / 02_plot_rmsd_rg_rmsf.py). Parse the
    .xvg by ROW POSITION, not the printed number: rows 1..chainA_length are
    chain A, the rest are chain B.
  - The reference PDB's residue numbering may be offset by a fixed amount
    from the simulation's own numbering (here, 4LNK residue N corresponds to
    simulation residue N-1 -- verify this for your own structure, it is NOT
    guaranteed to be 1).
  - A near-zero or even slightly negative Pearson r between B_sim and B_xray
    is a real possible outcome, not necessarily a bug -- crystallographic
    B-factors reflect both genuine flexibility AND crystal-packing/refinement
    effects that don't always track solution-phase MD dynamics. Confirmed in
    this project by checking that the correlation was equally weak in the
    pre-PBC-fix data, i.e. unrelated to the PBC bug that was fixed elsewhere.
"""
import sys
import json
import math
from collections import defaultdict


def parse_bfactors(pdb_path, resid_offset=-1):
    """CA-atom B-factors from a (possibly multi-chain/multi-copy) reference
    PDB, averaged per residue across all copies/chains present. resid_offset
    is ADDED to the PDB's own resid to map onto your simulation's numbering
    (default -1: reference residue N -> simulation residue N-1)."""
    bfac_by_resid = defaultdict(list)
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                resid = int(line[22:26])
                bfac = float(line[60:66])
                bfac_by_resid[resid].append(bfac)
    return {resid + resid_offset: sum(v) / len(v) for resid, v in bfac_by_resid.items()}


def parse_rmsf_by_chain(xvg_path, chainA_length=443):
    rmsf_values = []
    with open(xvg_path) as f:
        for line in f:
            if line.startswith(("#", "@")):
                continue
            parts = line.split()
            rmsf_values.append(float(parts[1]))
    chainA = {i + 1: v for i, v in enumerate(rmsf_values[0:chainA_length])}
    chainB = {i + 1: v for i, v in enumerate(rmsf_values[chainA_length:])}
    return chainA, chainB


def rmsf_to_bfactor(rmsf_nm):
    rmsf_A = rmsf_nm * 10.0
    return (8 * math.pi ** 2 / 3.0) * rmsf_A ** 2


def pearson(x, y):
    n = len(x)
    mx, my = sum(x) / n, sum(y) / n
    cov = sum((a - mx) * (b - my) for a, b in zip(x, y))
    vx = sum((a - mx) ** 2 for a in x)
    vy = sum((b - my) ** 2 for b in y)
    return cov / math.sqrt(vx * vy)


def compare(rmsf_xvg, ref_pdb, out_json, chainA_length=443, resid_offset=-1):
    bfac_avg = parse_bfactors(ref_pdb, resid_offset)
    chainA, chainB = parse_rmsf_by_chain(rmsf_xvg, chainA_length)

    out = {"chainA": [], "chainB": []}
    for label, chain_rmsf in (("chainA", chainA), ("chainB", chainB)):
        xs, ys = [], []
        for resid in sorted(chain_rmsf):
            if resid not in bfac_avg:
                continue
            rmsf_nm = chain_rmsf[resid]
            b_sim = rmsf_to_bfactor(rmsf_nm)
            b_xray = bfac_avg[resid]
            out[label].append({"resid": resid, "rmsf_nm": round(rmsf_nm, 4),
                                "b_sim": round(b_sim, 2), "b_xray": round(b_xray, 2)})
            xs.append(b_sim)
            ys.append(b_xray)
        r = pearson(xs, ys) if xs else float("nan")
        print(f"{label}: n={len(xs)} residues matched, Pearson r(B_sim, B_xray) = {r:.3f}")

    with open(out_json, "w") as f:
        json.dump(out, f)
    print("wrote", out_json)


if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: 05_rmsf_vs_bfactor.py <rmsf.xvg> <reference.pdb> <out.json> [chainA_length] [resid_offset]")
        print("Example (Run 1, as actually run):")
        print("  05_rmsf_vs_bfactor.py rmsf.xvg 4lnk.pdb rmsf_vs_4lnk.json 443 -1")
        sys.exit(1)
    rmsf_xvg, ref_pdb, out_json = sys.argv[1:4]
    chainA_length = int(sys.argv[4]) if len(sys.argv) > 4 else 443
    resid_offset = int(sys.argv[5]) if len(sys.argv) > 5 else -1
    compare(rmsf_xvg, ref_pdb, out_json, chainA_length, resid_offset)
