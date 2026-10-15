"""Report potential energy for structures already minimised, without re-running.

The minimisation log prints energy before and after. Where that log is lost,
this recomputes both from the saved files, using the same force field and
solvent model as minimise.py so the numbers are comparable to the ones
already reported.

    python3 energies.py gs_wt_AB_prot.pdb gs_wt_AB_min.pdb

Pass the pre-minimisation file first and the minimised file second. Any
number of files is accepted; each is reported on its own line.

The pre-minimisation file has no hydrogens and they are added here, exactly
as minimise.py does. The minimised file already has them, so addHydrogens is
skipped. Getting this wrong is the one way to produce a misleading number:
adding a second set of hydrogens to an already-protonated structure gives a
large spurious energy.
"""

import sys

from openmm import unit
from openmm.app import ForceField, Modeller, NoCutoff, PDBFile, Simulation
from openmm import LangevinIntegrator

FORCEFIELD = 'amber14-all.xml'
SOLVENT = 'implicit/obc2.xml'
PH = 7.0


def has_hydrogens(topology):
    return any(a.element is not None and a.element.symbol == 'H'
               for a in topology.atoms())


def energy(path):
    pdb = PDBFile(path)
    ff = ForceField(FORCEFIELD, SOLVENT)
    modeller = Modeller(pdb.topology, pdb.positions)

    if has_hydrogens(modeller.topology):
        added = False
    else:
        modeller.addHydrogens(ff, pH=PH)
        added = True

    system = ff.createSystem(modeller.topology, nonbondedMethod=NoCutoff,
                             constraints=None)
    integrator = LangevinIntegrator(300 * unit.kelvin, 1 / unit.picosecond,
                                    0.002 * unit.picoseconds)
    sim = Simulation(modeller.topology, system, integrator)
    sim.context.setPositions(modeller.positions)
    e = sim.context.getState(getEnergy=True).getPotentialEnergy()

    return (e.value_in_unit(unit.kilojoule_per_mole),
            modeller.topology.getNumAtoms(), added)


def main(paths):
    print(f"{'file':34}{'energy (kJ/mol)':>18}{'atoms':>9}  hydrogens")
    results = []
    for p in paths:
        try:
            e, n, added = energy(p)
            results.append(e)
            print(f"{p:34}{e:18.2f}{n:9d}  {'added here' if added else 'already present'}")
        except Exception as exc:
            print(f"{p:34}{'FAILED':>18}         {exc}")

    if len(results) == 2:
        print()
        print(f"  before : {results[0]:12.2f} kJ/mol")
        print(f"  after  : {results[1]:12.2f} kJ/mol")
        print(f"  change : {results[1] - results[0]:12.2f} kJ/mol")


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    main(sys.argv[1:])
