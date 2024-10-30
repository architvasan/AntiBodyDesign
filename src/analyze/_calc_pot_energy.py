from simtk.openmm.app import PDBFile
import MDAnalysis as mda
from simtk.openmm import app
from simtk.openmm import openmm as mm
from simtk import unit
import numpy as np
import tempfile

def main(inppdb, 
        inptraj,
        selA,
        selB):

    # Load the Universe
    u = mda.Universe(inppdb, inptraj)
    
    # Define force field and platform
    forcefield = app.ForceField('amber14-all.xml', 'tip3p.xml')
    platform = mm.Platform.getPlatformByName('CPU')
    
    # Define the pair of residues for interaction energy calculation
    res1 = u.select_atoms(selA)
    res2 = u.select_atoms(selB)
    
    # Set up the integrator
    integrator = mm.VerletIntegrator(0.001*unit.picoseconds)
    
    # Iterate over the trajectory
    for ts in u.trajectory:
        # Select the positions of the two residues
        combined_sel = u.select_atoms(f'({selA}) or ({selB})')
        
        # Write the selected atoms to a temporary PDB file
        with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as tmpfile:
            combined_sel.write(tmpfile.name)
            from pdbfixer import PDBFixer
            # Load the PDB file
            fixer = PDBFixer(filename=tmpfile.name)
            
            # Add missing residues, chains, and atoms
            fixer.findMissingResidues()
            fixer.findNonstandardResidues()
            fixer.replaceNonstandardResidues()
            fixer.findMissingAtoms()
            fixer.addMissingAtoms()
            fixer.addMissingHydrogens()
            # Add missing termini (important for resolving C-terminal issues)
            #fixer.missingTerminals()

            # Save the fixed PDB file
            with open(tmpfile.name, 'w') as f:
                PDBFile.writeFile(fixer.topology, fixer.positions, f)

            pdbfile = app.PDBFile(tmpfile.name)
        
        # Create OpenMM system for this subsystem
        system = forcefield.createSystem(pdbfile.topology, nonbondedMethod=app.NoCutoff, constraints=None)
    
        # Create a simulation context
        simulation = app.Simulation(pdbfile.topology, system, integrator, platform)
        simulation.context.setPositions(combined_sel.positions)
    
        # Get potential energy of the combined system
        state = simulation.context.getState(getEnergy=True)
        energy_total = state.getPotentialEnergy()
    
        # Calculate interaction energy by excluding the interactions within each residue
        # Get energy of res1 alone
        simulation.context.setPositions(res1.positions)
        state_res1 = simulation.context.getState(getEnergy=True)
        energy_res1 = state_res1.getPotentialEnergy()
    
        # Get energy of res2 alone
        simulation.context.setPositions(res2.positions)
        state_res2 = simulation.context.getState(getEnergy=True)
        energy_res2 = state_res2.getPotentialEnergy()
    
        # Interaction energy = Total energy - (Energy of res1 alone + Energy of res2 alone)
        interaction_energy = energy_total - (energy_res1 + energy_res2)
    
        print(f'Frame {ts.frame}: Interaction Energy = {interaction_energy}')

if __name__ == "__main__":

    def list_of_strings(arg):
        return arg.split(',')

    import argparse
    parser = argparse.ArgumentParser()

    parser.add_argument('-p',
                        '--inputpdb',
                        type=str,
                        help='input pdb with protein')

#    parser.add_argument('-t',
#                        '--inputtraj',
#                        required=False,
#                        type=str,
#                        default='none',
#                        help='single trajectory file')

    parser.add_argument('-T',
                        '--inputtraj_list',
                        required=False,
                        type=list_of_strings,
                        default=['none', 'none'],
                        help='trajectory file list (dcd format)')

    parser.add_argument('-sA',
                        '--selA',
                        type=str,
                        help='phrase for seelction A (in mdanalysis language)')

    parser.add_argument('-sB',
                        '--selB',
                        type=str,
                        help='phrase for seelction A (in mdanalysis language)')

    args = parser.parse_args()
    main(args.inputpdb, args.inputtraj_list, args.selA, args.selB)
