from pdbfixer import PDBFixer
from openmm.app import PDBFile

def main(inputpdb, 
         outputpdb):
    # Load the PDB file
    fixer = PDBFixer(filename=inputpdb)
    
    # Iterate over the residues and change HIE to HSE
    for residue in fixer.topology.residues():
        if residue.name == 'HIE':
            residue.name = 'HSE'
    
    # Save the modified PDB file
    with open(outputpdb, 'w') as output_file:
        PDBFile.writeFile(fixer.topology, fixer.positions, output_file)

if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser()

    parser.add_argument('-p',
                        '--inputpdb',
                        type=str,
                        help='input pdb with protein')

    parser.add_argument('-o',
                         '--outputpdb',
                         type=str,
                         help='output cleaned pdb')

    args = parser.parse_args()
    main(args.inputpdb, args.outputpdb)
