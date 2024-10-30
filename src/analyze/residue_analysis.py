import MDAnalysis as mda

def residuecount(inputpdb):
    # Load your PDB file
    u = mda.Universe(inputpdb)
    # Select atoms based on the criteria
    selected_residues = u.select_atoms('segid A').residues
    print(selected_residues)
    # Define residue property categories
    positively_charged_residues = {'ARG', 'LYS'}
    negatively_charged_residues = {'ASP', 'GLU'}
    polar_residues = {'ASN', 'GLN', 'SER', 'THR', 'TYR'}
    nonpolar_residues = {'ALA', 'CYS', 'GLY', 'ILE', 'LEU', 'MET', 'PHE', 'PRO', 'TRP', 'VAL'}
    
    # Initialize counters
    pos_charged_count = 0
    neg_charged_count = 0
    polar_count = 0
    nonpolar_count = 0
    
    # Loop through residues in the PDB file
    for residue in selected_residues:
        res_name = residue.resname
        if res_name in positively_charged_residues:
            pos_charged_count += 1
        elif res_name in negatively_charged_residues:
            neg_charged_count += 1
        elif res_name in polar_residues:
            polar_count += 1
        elif res_name in nonpolar_residues:
            nonpolar_count += 1
    
    # Print the results
    #print(f'Positively charged residues: {pos_charged_count}')
    #print(f'Negatively charged residues: {neg_charged_count}')
    #print(f'Polar residues: {polar_count}')
    #print(f'Nonpolar residues: {nonpolar_count}')
    return pos_charged_count, neg_charged_count, polar_count, nonpolar_count

if __name__ == "__main__":

    def list_of_strings(arg):
        return arg.split(',')

    import argparse
    parser = argparse.ArgumentParser()

    parser.add_argument('-p',
                        '--inputpdb',
                        type=str,
                        help='input pdb with protein')

    parser.add_argument('-O',
                        '--outdir',
                        type=str,
                        help='directory to output data')

    parser.add_argument('-o',
                        '--outfile',
                        type=str,
                        help='file to output')


    args = parser.parse_args()

    pos_charged_count, neg_charged_count, polar_count, nonpolar_count = residuecount(args.inputpdb)

        # Write the results to the file
    with open(output_file, 'w') as file:
        file.write(f'Pos,Neg,Polar,NonPolar\n')
        file.write(f'{pos_charged_count},{neg_charged_count},{polar_count},{nonpolar_count}')
