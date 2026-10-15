"""Folding ddG for H195Y by replicated local cartesian repack-and-minimise.

Protocol: ref2015_cart, shell of residues within 8 A of the mutation site,
repack + cartesian minimise, n replicates, wild type and mutant from the same
starting pose under an identical protocol.

Sign convention here: ddG = mutant - wild type.
    positive  -> mutation is DESTABILISING
    negative  -> mutation is STABILISING
Divide REU by ~2.9 for an approximate kcal/mol.

    python3 ddg_h195y.py gs_wt_AB_prot.pdb A 195 TYR

Run on the A-B DIMER, not a monomer: the GS active site draws four residues
from the partner chain, and position 195 sits at the feedback site. A monomer
gives half a pocket.

Ligands and metals must be stripped first (PyRosetta has no params for ADP,
glutamate or Mg without supplying them). The comparison is therefore
differential and valid; absolute energies are not comparable to any other
structure prepared differently.
"""

import sys
import statistics as stats

from pyrosetta import init, pose_from_pdb, create_score_function
from pyrosetta.rosetta.core.kinematics import MoveMap
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.pack.task.operation import (
    InitializeFromCommandline, RestrictToRepacking, PreventRepacking,
)
from pyrosetta.rosetta.core.select.residue_selector import (
    ResidueIndexSelector, NeighborhoodResidueSelector,
)
from pyrosetta.rosetta.protocols.minimization_packing import (
    PackRotamersMover, MinMover,
)
from pyrosetta.rosetta.protocols.simple_moves import MutateResidue

SHELL_RADIUS = 8.0
N_REPLICATES = 5
MIN_TOLERANCE = 1e-4
MAX_MIN_ITER = 2000
REU_PER_KCAL = 2.9


def shell_indices(pose, target_idx, radius=SHELL_RADIUS):
    """Pose indices within `radius` of the target residue, target included."""
    sel = NeighborhoodResidueSelector(
        ResidueIndexSelector(str(target_idx)), radius, True
    )
    mask = sel.apply(pose)
    return [i for i in range(1, pose.total_residue() + 1) if mask[i]]


def build_movemap(pose, idx_list):
    mm = MoveMap()
    mm.set_bb(False)
    mm.set_chi(False)
    for i in idx_list:
        mm.set_bb(i, True)
        mm.set_chi(i, True)
    return mm


def relax_shell(pose, sf, idx_list):
    """Repack then cartesian-minimise the shell. Modifies pose in place."""
    tf = TaskFactory()
    tf.push_back(InitializeFromCommandline())
    tf.push_back(RestrictToRepacking())

    keep_fixed = PreventRepacking()
    inside = set(idx_list)
    for i in range(1, pose.total_residue() + 1):
        if i not in inside:
            keep_fixed.include_residue(i)
    tf.push_back(keep_fixed)

    PackRotamersMover(sf, tf.create_task_and_apply_taskoperations(pose)).apply(pose)

    mm = build_movemap(pose, idx_list)
    mn = MinMover(mm, sf, 'lbfgs_armijo_nonmonotone', MIN_TOLERANCE, True)
    mn.cartesian(True)
    mn.max_iter(MAX_MIN_ITER)
    mn.apply(pose)
    return sf(pose)


def one_replicate(pdb, sf, chain, resnum, mutant_aa):
    """Returns (wt_score, mut_score) from the same starting coordinates."""
    base = pose_from_pdb(pdb)
    idx = base.pdb_info().pdb2pose(chain, int(resnum))
    if idx == 0:
        raise SystemExit(f'residue {chain}/{resnum} not found in {pdb}')

    native = base.residue(idx).name3()
    shell = shell_indices(base, idx)

    wt = base.clone()
    wt_score = relax_shell(wt, sf, shell)

    mut = base.clone()
    MutateResidue(idx, mutant_aa).apply(mut)
    mut_score = relax_shell(mut, sf, shell)

    return wt_score, mut_score, native, idx, len(shell)


def main():
    if len(sys.argv) != 5:
        raise SystemExit(__doc__)
    pdb, chain, resnum, mutant_aa = sys.argv[1:5]

    init('-mute all -ignore_unrecognized_res -ignore_zero_occupancy false '
         '-ex1 -ex2 -use_input_sc -flip_HNQ')
    sf = create_score_function('ref2015_cart')

    wt_scores, mut_scores = [], []
    native = idx = n_shell = None

    for r in range(1, N_REPLICATES + 1):
        w, m, native, idx, n_shell = one_replicate(pdb, sf, chain, resnum, mutant_aa)
        wt_scores.append(w)
        mut_scores.append(m)
        print(f'  replicate {r}:  WT {w:10.2f}   mut {m:10.2f}   ddG {m-w:+8.2f} REU')

    ddg = [m - w for w, m in zip(wt_scores, mut_scores)]
    mean = stats.mean(ddg)
    sd = stats.stdev(ddg) if len(ddg) > 1 else 0.0

    print(f'\n=== {native}{resnum}{mutant_aa}, chain {chain} ===')
    print(f'  pose index         {idx}')
    print(f'  shell residues     {n_shell} within {SHELL_RADIUS:.0f} A')
    print(f'  replicates         {N_REPLICATES}')
    print(f'  WT   mean +- sd    {stats.mean(wt_scores):10.2f} +- '
          f'{stats.stdev(wt_scores) if len(wt_scores)>1 else 0:.2f} REU')
    print(f'  mut  mean +- sd    {stats.mean(mut_scores):10.2f} +- '
          f'{stats.stdev(mut_scores) if len(mut_scores)>1 else 0:.2f} REU')
    print(f'  ddG  mean +- sd    {mean:+10.2f} +- {sd:.2f} REU')
    print(f'  ddG  approx        {mean/REU_PER_KCAL:+10.2f} kcal/mol')
    print(f'\n  {"DESTABILISING" if mean > 0 else "STABILISING"} by this protocol.')

    if sd > abs(mean):
        print('\n  WARNING: spread exceeds the mean. Raise N_REPLICATES or report'
              '\n  the result as indistinguishable from zero.')
    if abs(mean) < 1.0:
        print('\n  NOTE: |ddG| below 1 REU (~0.34 kcal/mol) is within the noise of'
              '\n  this protocol. Report as no detectable effect rather than as a'
              '\n  small stabilisation or destabilisation.')


if __name__ == '__main__':
    main()
