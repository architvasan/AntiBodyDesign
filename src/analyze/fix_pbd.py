from pdbfixer import PDBFixer
from openmm.app import PDBFile
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('-p',
                    '--inputpdb',
                    type=str,
                    help='input pdb with protein')

parser.add_argument('-o',
                    '--outputpdb',
                    type=str,
                    help='input pdb with protein')





args = parser.parse_args()

# Load your PDB file with PDBFixer
fixer = PDBFixer(filename=args.inputpdb)

# Find missing residues and atoms
fixer.findMissingResidues()
fixer.findMissingAtoms()
fixer.addMissingAtoms()

# Add missing hydrogens (assuming pH 7.0)
fixer.addMissingHydrogens(pH=7.0)

# Save the fixed PDB file
with open(args.outputpdb, 'w') as output:
    PDBFile.writeFile(fixer.topology, fixer.positions, output)

