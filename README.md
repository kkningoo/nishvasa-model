# Nishvasa — model code

Code behind the Model page for iGEM VIT Vellore 2026. An enzymatic ammonia
scrubber for poultry sheds: glutamine synthetase H195Y and polyphosphate
kinase 2 co-immobilised in chitosan beads inside a 3D-printed PETG honeycomb.

`anchormap`, the immobilisation and interface analysis as an installable
package, is at https://github.com/kkningoo/anchormap. It is a cleaner
re-implementation of much of `gate2_structural/`. These scripts are what
actually produced the numbers on the page; the package is what another team
should use.

## Layout

The directories follow the four gates of the Model page.

```
parameters.yaml                every value with its provenance class

gate1_kinetics/                which enzymes, in which variants
  model1_gs_feedback.py        Cheng-Prusoff, four-term rate, wild type vs
                               H195Y, and the two Gate 1 figures
  variant_selection.py         all six variants, ratios and scatter

gate2_structural/              can they be made into one component
  prep/
    minimise.py                restrained minimisation, OpenMM + AMBER14/OBC2
    energies.py                single-point energies from saved structures
  accessibility/
    gs_surface.py              GS lysines, both oligomeric states, + terminals
    ppk2_sasa.py               PPK2 per-lysine accessibility
    ppk2_sasa2.py              PPK2 interface areas by FreeSASA
    liabilities.py             deamidation, isomerisation, oxidation sites
  coupling/
    multipoint.py              PPK2 multipoint reachability
    multipoint2.py             GS multipoint reachability
    ringspan.py                ring-stratified span sweep, both floors
  interfaces/
    interface.py               GS pairwise interface energetics
    interface2.py              GS C-D variant
    interface_ppk2.py          PPK2, backbone fixed
    interface_ppk2_cd.py       PPK2 C-D, backbone free -- NOT comparable
    interfaces.py              FreeSASA buried areas, all pairs
    ddg_h195y.py               cartesian ddG, H195Y
  electrostatics/
    pka_summary.py             PROPKA output parsing
    plot_pka.py                stability curves, metals retained
    plot_pka_no_mg.py          stability curves, metals stripped
  figures/
    interface_figure.py        energy density and the stitch-count reversal
  md/                          [to be added]

gate3_coupled/                 does the coupled reaction sustain itself
  model2_atp_regeneration.py   2A coating liquid, PPK2 necessity
                               2B reservoir-coupled 7-day run

gate4_device/                  can it be delivered at a useful rate
  model3_bead_effectiveness.py Thiele modulus, three shapes, diffusivity audit
  model4_pegda_coating.py      coating thickness sweep in the real cell
  model5_scaffold_geometry.py  hexagonal cell inradius, VRE and design score
  model6_cartridge_removal.py  plug-flow integration along the channel
  model7_batch_scrubbing.py    batch time and farm sizing
  model8_continuous_operation.py  continuous mode, rejected at 7-70 units

gate5_sensitivity/             [to be added]

data/                          structures and parsed outputs [not committed]
```

## Running

Most scripts take their inputs at the top of the file and are run directly:

```bash
python3 gate1_kinetics/model1_gs_feedback.py
python3 gate4_device/model6_cartridge_removal.py
```

`model6` and `model7` integrate along the channel at every concentration and
take a few minutes. Everything else is seconds.

Dependencies: numpy, scipy, matplotlib for Gates 1, 3 and 4; biopython and
freesasa for the Gate 2 accessibility work. The interface scripts need
PyRosetta, which is not pip-installable and requires an academic licence.
The minimisation needs OpenMM. `requirements.txt` lists both groups.

## Reproducing the page

| Script | Produces |
| --- | --- |
| `model1_gs_feedback.py` | Ki 19.2 vs 0.37 mM, 52-fold; H195Y clears 99.9 % in 60 min, wild type stalls at 77.2 % |
| `model2_atp_regeneration.py` | ATP collapse without PPK2; phosphate crosses 12.5 mM 2.4x before polyphosphate is spent |
| `model3_bead_effectiveness.py` | eta_bead = 0.9995 at 10 um; the withdrawn tuned constant implied a diffusivity 40-100x too low |
| `model4_pegda_coating.py` | eta_coat = 0.985 at 24 um; J peaks at 256 um, which is not why 24 um was chosen |
| `model5_scaffold_geometry.py` | r = 1.5 mm at 99.2 % of the VRE peak; the circumscribed-circle error was a 21 % overestimate |
| `model6_cartridge_removal.py` | 9,244 cells, 26.5 m2, 713 umol/min at 50 ppm, 13.6 Pa |
| `model7_batch_scrubbing.py` | 6.0 / 7.3 / 9.6 h for the three farm sizes |
| `model8_continuous_operation.py` | 7-70 units continuous against 1-6 batch |

## Three corrections carried in this copy

These scripts were corrected after producing the numbers they produced. The
corrections are in the code; the reasoning is on the Model page.

**`interface_ppk2.py` -- `set_bb(False)` restored.** The tabled A-C and A-B
values were produced with the backbone fixed. The file was later edited in
place when the C-D retry was run backbone-free, so the archived copy no longer
matched the protocol behind its own results. `interface_ppk2_cd.py`
deliberately differs and its result is reported separately, not ranked.

**`multipoint.py`, `multipoint2.py`, `ringspan.py` -- span floor 10 A to 5 A.**
A bifunctional crosslinker is limited by how far it can reach, not by how
close two amines may be, and genipin self-polymerises into short dimers and
trimers. The 10 A floor excluded four of six real cross-ring contacts at
8.4-8.5 A. All three now use 5-16 A, and `ringspan.py` reports both floors so
the dependence is visible.

**`liabilities.py` -- Tien theoretical maxima.** The original used MET 203,
SER 130 and THR 142 from a different table, inflating methionine RSA by about
10 % and threonine by about 21 %. No liability crosses the cutoff under the
correction; only the numbers move.

## What the parameters file is for

`parameters.yaml` carries every value on the Model page with one of four
provenance classes: literature, prediction, assumption, or measurement. The
measurement list is empty, deliberately, and the file says so -- no wet-lab
measurement was available this season and every rate inherits that.

Withdrawn values are kept with a note rather than deleted. The effectiveness
factor that was 0.70 and is now 1.00, the removal rate that was 1022 and is
now 715, the pressure drop that was 0.91 Pa and is now 13.6, the span floor
that was 10 A and is now 5 A: each is a correction worth being able to trace.

Three values are quoted at different precision on the page and in the Gate 4
working -- Ki 19.2 against 19.17, airflow 180 against 181, pressure drop 14 Pa
against 13.6. The file carries the precise value and notes the rounded one.

Two entries carry an unverified source and say so: Km,PolyP for PPK2, and the
ammonia Michaelis constant whose Nocek attribution is withdrawn.

## Licence

MIT. See `LICENSE`.
