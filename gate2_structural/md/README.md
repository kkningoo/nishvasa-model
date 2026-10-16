# Molecular dynamics

Three runs on the GS A-B dimer, reported in Model page section 2.6 and
supplementary S5b.

| Run | System | Deliverable |
| --- | --- | --- |
| 1 | Wild type holo, metals and ligands retained | RMSF against crystallographic B-factors |
| 2 | Wild type and H195Y, inhibited state | Contact lifetimes at the glutamine site |
| 3 | Wild type, magnesium depleted | Interface contact retention |

## Order

Run `01` first on any trajectory, before trusting its RMSD or RMSF. The
remaining scripts assume a PBC-corrected trajectory.

```
01_pbc_fix_and_rmsd_rmsf_gyrate.sh   two-pass trjconv, then RMSD, RMSF, Rg, H-bonds
02_plot_rmsd_rg_rmsf.py              plots the .xvg output, RMSF split by chain
03_contact_lifetime_analysis.py      template, with the gotchas documented
03b_..._WORKED_EXAMPLE_wt_inhib.py   the filled-in version that produced Run 2
04_interface_contact_stability.py    native inter-chain contact retention (Run 3)
05_rmsf_vs_bfactor.py                RMSF against crystallographic B-factors (Run 1)
```

Requires `gmx` and MDAnalysis.

## The PBC bug, and what it means for earlier figures

Every one of these runs' original trajectories had the same defect: a naive
`trjconv -pbc nojump` on a multi-chain system unwraps each molecule
independently without first clustering the assembly, so over a long
trajectory the chains drift into different periodic images relative to each
other. Backbone RMSD blows up to tens of nanometres and per-residue RMSF
produces 5-12 nm spikes at the chain boundary, while the physical structure
is perfectly stable -- the radius of gyration stays flat throughout, which is
the tell.

`01` fixes it with a two-pass `trjconv`: `-pbc cluster` on a non-water
selection, then `-pbc nojump` on the result. Clustering alone is not enough,
because per-frame clustering is not temporally stable and can flip between
equally valid periodic images, which makes RMSD oscillate.

A second defect compounds it. `gmx rms` and `gmx rmsf` take their fit
reference from the coordinates embedded in the `-s` tpr, not from the
trajectory. A tpr built from a mid-production restart is very likely in a
different periodic image from the corrected trajectory, which reintroduces a
large constant offset. `01` extracts frame 0 of the corrected trajectory and
builds a throwaway tpr from it, so the reference is self-consistent.

**This is why the two 100 ns RMSD traces reported on the page start from
different values.** Both are affected; they should be regenerated with `01`
before the plots are used for anything quantitative.

## Three things these scripts settle

**The residue-numbering offset.** The page uses Asp53, Arg62 and Glu65; these
scripts use Asp52, Arg61 and Glu64. The simulation numbers residues
continuously across both chains, 1-886 with no restart, and chain B's
equivalent of chain A's residue X is X+443, not X+444. The simulation's own
numbering runs one below the crystallographic numbering. Both sets refer to
the same residues.

**Which magnesium is which.** There are four Mg ions, two per chain. The
catalytic site, n2, is identified by measuring which sits about 2 A from the
bound glutamine OE1 -- resid 890 for chain A, 892 for chain B. The other two,
891 and 893, are the n1 site near the phosphate. Not identifiable by file
order.

**Why the B-factor correlation is weak.** A near-zero or slightly negative
Pearson r between simulated and crystallographic B-factors is a real outcome,
not a bug. Crystallographic B-factors carry lattice packing and refinement
effects alongside genuine flexibility. Confirmed here by checking that the
correlation was equally weak before the PBC fix, so it is unrelated to that
defect.

## What is not here

The trajectories themselves and the .tpr files. They are too large to commit
and are not reproducible from this repository alone. The paths at the top of
`03b` and `04` are the machine they were run on; edit before reuse.
