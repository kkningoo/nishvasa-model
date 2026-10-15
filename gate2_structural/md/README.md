# Molecular dynamics

To be added. Three runs, reported in section 2.6 and supplementary S5b.

| Run | System | Length |
| --- | --- | --- |
| Holo | GS wild type, A-B dimer, metals and ligands retained | 100 ns, 308 K |
| Mg-depleted | GS wild type, A-B dimer, magnesium removed | 100 ns, 308 K |
| H195Y inhibited | GS H195Y, 4LNF-derived, inhibited state | 10 ns |

Needed here: the run configuration for each, and the analysis scripts
behind the RMSD, radius of gyration, RMSF, native-contact and
B-factor-comparison plots.

Two things flagged as unresolved on the page that these scripts should
make checkable:

- the two 100 ns RMSD traces use different reference structures, so their
  absolute values are not comparable as recorded
- the H195Y run reports Asp52, Arg61 and Glu64 where the rest of this work
  uses Asp53, Arg62 and Glu65, a one-residue offset most likely from
  PDBFixer renumbering
