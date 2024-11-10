from pdbfixer import PDBFixer
from openmm import *
from openmm.app import *
from openmm.unit import *
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

if __name__ == "__main__": 
    import argparse
    parser = argparse.ArgumentParser()
    
    parser.add_argument('-p',
                        '--pdbinp',
                        type=str,
                        help='input pdb')

    parser.add_argument('-op',
                        '--outpdb',
                        type=str,
                        help='input pdb')

    args = parser.parse_args()

    pdbfixit(args.pdbinp, args.outpdb)