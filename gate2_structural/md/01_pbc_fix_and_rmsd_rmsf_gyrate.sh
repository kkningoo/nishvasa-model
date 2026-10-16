#!/bin/bash
# PBC-correct a GROMACS trajectory and compute RMSD, RMSF, radius of gyration,
# and H-bond count on the corrected trajectory.
#
# Background: a naive `gmx trjconv -pbc nojump` on a multi-chain system (e.g. a
# homodimer + separate ligand molecules) unwraps each molecule's own PBC jumps
# independently, without keeping the whole assembly together. Over a long
# trajectory this lets chains drift into different periodic images relative to
# each other, which blows up backbone RMSD (seen as high as 18-35 nm) and
# produces impossible per-residue RMSF spikes (5-12 nm) at whichever residue
# happens to sit at the chain-boundary in the fit group's internal numbering --
# even though the real, physical structure is perfectly stable (radius of
# gyration stays flat the whole time).
#
# The fix used here is a two-pass trjconv:
#   1. `-pbc cluster` on a protein+ligands+ions selection (NOT System -- do not
#      include bulk water, clustering tens of thousands of water molecules is
#      needlessly expensive and can take hours) to join everything into one
#      periodic image per frame. Output group must still be "System" so atom
#      count/indexing matches the original topology for later analysis steps.
#   2. `-pbc nojump` on top of that, because per-frame clustering is NOT
#      temporally stable by itself -- it can flip between equally-valid
#      periodic images frame-to-frame, which (counterintuitively) makes RMSD
#      *worse* (oscillating) if you skip this second pass.
#
# A second, separate gotcha: `gmx rms`/`gmx rmsf` take their fit REFERENCE
# structure from whatever coordinates are embedded in the `-s` tpr file, not
# from the trajectory itself. If that tpr was built from an arbitrary
# mid-production restart point (e.g. a segment's own grompp invocation), its
# embedded structure is very likely in a DIFFERENT periodic image than your
# freshly-corrected trajectory, which reintroduces a large (if now constant,
# non-exploding) RMSD/RMSF offset. The fix: extract frame 0 of your corrected
# trajectory and build a throwaway tpr from it (via grompp with a minimal mdp)
# so the reference is self-consistent with the trajectory you're analyzing.
#
# Usage:
#   01_pbc_fix_and_rmsd_rmsf_gyrate.sh <workdir> <tpr> <xtc> <topol_top> <minim_mdp> <tag>
#
# Requires in <workdir>: the tpr/xtc, topol.top, and whatever .itp files
# topol.top #includes (posre_*.itp etc), plus a minimal minim.mdp (just needs
# to be grompp-able; its integrator/nsteps are irrelevant since nothing is run).
#
# Outputs (in <workdir>): <tag>_fixed_nojump.xtc, rmsd_<tag>.xvg,
# rmsf_<tag>.xvg, gyrate_<tag>.xvg, hbnum_<tag>.xvg, ref0_<tag>.gro/.tpr

set -e
GMX=${GMX:-gmx}

WORKDIR="$1"; TPR="$2"; XTC="$3"; TOP="$4"; MINIM_MDP="$5"; TAG="$6"
if [ -z "$TAG" ]; then
  echo "Usage: $0 <workdir> <tpr> <xtc> <topol_top> <minim_mdp> <tag>"
  exit 1
fi
cd "$WORKDIR"

echo "=== [$TAG] Step 1: cluster (non-Water group) then nojump ==="
printf "non-Water\nSystem\n" | "$GMX" trjconv -s "$TPR" -f "$XTC" -pbc cluster -o "${TAG}_clustered.xtc"
printf "System\n" | "$GMX" trjconv -s "$TPR" -f "${TAG}_clustered.xtc" -pbc nojump -o "${TAG}_fixed_nojump.xtc"
rm -f "${TAG}_clustered.xtc"

echo "=== [$TAG] Step 2: self-consistent reference (frame 0 of corrected trajectory) ==="
printf "System\n" | "$GMX" trjconv -s "$TPR" -f "${TAG}_fixed_nojump.xtc" -dump 0 -o "ref0_${TAG}.gro"
"$GMX" grompp -f "$MINIM_MDP" -c "ref0_${TAG}.gro" -p "$TOP" -o "ref0_${TAG}.tpr" -maxwarn 20

echo "=== [$TAG] Step 3: RMSD (backbone, lsq fit) ==="
printf "Backbone\nBackbone\n" | "$GMX" rms -s "ref0_${TAG}.tpr" -f "${TAG}_fixed_nojump.xtc" -o "rmsd_${TAG}.xvg" -tu ns

echo "=== [$TAG] Step 4: RMSF (per residue, fit to backbone) ==="
printf "Backbone\nProtein\n" | "$GMX" rmsf -s "ref0_${TAG}.tpr" -f "${TAG}_fixed_nojump.xtc" -o "rmsf_${TAG}.xvg" -res -fit

echo "=== [$TAG] Step 5: radius of gyration ==="
printf "Protein\n" | "$GMX" gyrate -s "ref0_${TAG}.tpr" -f "${TAG}_fixed_nojump.xtc" -o "gyrate_${TAG}.xvg"

echo "=== [$TAG] Step 6: H-bond count ==="
printf "Protein\nProtein\n" | "$GMX" hbond -s "ref0_${TAG}.tpr" -f "${TAG}_fixed_nojump.xtc" -num "hbnum_${TAG}.xvg"

echo "=== [$TAG] DONE -- sanity check these before trusting them: ==="
echo "--- RMSD head/tail (should be a smooth, bounded curve -- a few tenths to a few nm over 100ns, NOT tens of nm) ---"
grep -v "^[#@]" "rmsd_${TAG}.xvg" | head -3
grep -v "^[#@]" "rmsd_${TAG}.xvg" | tail -3
echo "--- RMSF top5 (watch for an isolated spike at the chain-A/chain-B boundary row -- if present after this fix, it's likely real terminal/ligand flexibility, not PBC) ---"
grep -v "^[#@]" "rmsf_${TAG}.xvg" | sort -k2 -n -r | head -5
echo "--- Rg head/tail (should stay roughly flat if the complex is physically stable) ---"
grep -v "^[#@]" "gyrate_${TAG}.xvg" | head -3
grep -v "^[#@]" "gyrate_${TAG}.xvg" | tail -3
