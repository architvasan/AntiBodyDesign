from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from tqdm import tqdm
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('-p',
                    '--inputpdb',
                    type=str,
                    help='input pdb with protein')

args = parser.parse_args()

# Load the PDB file
pdb = PDBFile(args.inputpdb)

# Load force field parameters
forcefield = ForceField('amber14-all.xml', 'amber14/tip3p.xml')

# Create a system with the force field
system = forcefield.createSystem(pdb.topology, nonbondedMethod=NoCutoff, constraints=HBonds)

# Load force field parameters
forcefield = ForceField('amber14-all.xml', 'amber14/tip3p.xml')

# Create a system with the force field
system = forcefield.createSystem(pdb.topology, nonbondedMethod=NoCutoff, constraints=HBonds)

# Create a context to perform energy calculation
integrator = VerletIntegrator(0.001*picoseconds)
context = Context(system, integrator)
context.setPositions(pdb.positions)

# Calculate the total nonbonded interaction energy
state = context.getState(getEnergy=True)
interaction_energy = state.getPotentialEnergy()

print(f'Total Interaction Energy: {interaction_energy}')
