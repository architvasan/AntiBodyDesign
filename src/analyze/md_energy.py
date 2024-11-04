from simtk.openmm import app, openmm
from simtk import unit
from openmm.app import PDBFile


def pdbfixit(inpdb, outpdb):
    fixer = PDBFixer(filename=inpdb)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens()
    
    # Save the fixed structure to a temporary file
    fixed_pdb_file = outpdb
    with open(fixed_pdb_file, 'w') as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f)


def calculate_md_energy(pdb_file):
    # Load the PDB structure
    pdb = PDBFile(pdb_file)
    
    # Define the force field
    forcefield = app.ForceField('amber14-all.xml', 'tip3p.xml')
    
    # Create the system from the topology
    system = forcefield.createSystem(
        pdb.topology,
        nonbondedMethod=app.NoCutoff,  # No periodic boundaries for single structure energy calc
        constraints=app.HBonds
    )
    
    # Set up an integrator (does not affect energy calculation here)
    integrator = openmm.VerletIntegrator(0.001 * unit.picoseconds)
    
    # Create a simulation context
    platform = openmm.Platform.getPlatformByName('CUDA')
    simulation = app.Simulation(pdb.topology, system, integrator, platform)
    
    # Set the particle positions
    simulation.context.setPositions(pdb.positions)
    
    # Minimize the energy
    simulation.minimizeEnergy()
    
    # Get the state and print the potential energy
    state = simulation.context.getState(getEnergy=True)
    potential_energy = state.getPotentialEnergy()
    
    print(f"Potential Energy: {potential_energy.value_in_unit(unit.kilocalories_per_mole)} kcal/mol")
    return potential_energy.value_in_unit(unit.kilocalories_per_mole)

if __name__ == "__main__": 
    import argparse
    parser = argparse.ArgumentParser()
    
    parser.add_argument('-p',
                        '--pdbinp',
                        type=str,
                        help='input pdb')

    args = parser.parse_args()
    # Usage
    calculate_md_energy(args.pdbinp)  # Replace 'your_structure.pdb' with the path to your PDB file

