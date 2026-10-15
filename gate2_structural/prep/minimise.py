from openmm.app import *
from openmm import *
from openmm.unit import *

pdb = PDBFile('4lnf_mut_prot.pdb')   # ADP and Mg stripped
ff = ForceField('amber14-all.xml', 'implicit/obc2.xml')

modeller = Modeller(pdb.topology, pdb.positions)
modeller.addHydrogens(ff, pH=7.0)

system = ff.createSystem(modeller.topology,
                         nonbondedMethod=NoCutoff,
                         constraints=None)

# restrain heavy atoms
force = CustomExternalForce('k*periodicdistance(x,y,z,x0,y0,z0)^2')
force.addGlobalParameter('k', 1000.0*kilojoules_per_mole/nanometer**2)
for p in ['x0','y0','z0']:
    force.addPerParticleParameter(p)
for atom in modeller.topology.atoms():
    if atom.element.symbol != 'H':
        pos = modeller.positions[atom.index]
        force.addParticle(atom.index, [pos.x, pos.y, pos.z])
system.addForce(force)

integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
sim = Simulation(modeller.topology, system, integrator)
sim.context.setPositions(modeller.positions)

print('E before:', sim.context.getState(getEnergy=True).getPotentialEnergy())
sim.minimizeEnergy(maxIterations=2000)
print('E after: ', sim.context.getState(getEnergy=True).getPotentialEnergy())

state = sim.context.getState(getPositions=True)
PDBFile.writeFile(modeller.topology, state.getPositions(),
                  open('4lnf_mut_min.pdb','w'))
