import pymol2

def _truncate_protein_to_nearby_residues(pdb_file, selection, output_file, distance_cutoff=15.0):
    with pymol2.PyMOL() as pymol:
        pymol.cmd.reinitialize()
        
        # Load the PDB file
        pymol.cmd.load(pdb_file, "protein")

        # Define the selection and specify residues within the cutoff distance
        pymol.cmd.select("close_residues", f"byres (br. ({selection}) within {distance_cutoff})")
        
        # Save the truncated selection to a new file
        pymol.cmd.save(output_file, "close_residues")
        
        print(f"Truncated protein saved to {output_file}")

def _truncate_pdb(pdb_file, selection, output_file, distance_cutoff=15.0):
    with pymol2.PyMOL() as pymol:
        pymol.cmd.reinitialize()
        
        # Load the PDB file
        pymol.cmd.load(pdb_file, "protein")
        
        # Select residues within the distance cutoff
        pymol.cmd.select("close_residues", f"byres (br. ({selection}) within {distance_cutoff})")
        
        # Extract and save the selection to a temporary object
        pymol.cmd.create("truncated_protein", "close_residues")
        
        # Get original residue IDs and store them for mapping
        original_residues = []
        pymol.cmd.iterate("truncated_protein", "original_residues.append(resi)", space={'original_residues': original_residues})
        
        # Renumber residues contiguously starting from 1
        pymol.cmd.alter("truncated_protein", "resi=str(int(resi))")  # Normalize residue IDs to integers
        pymol.cmd.sort("truncated_protein")  # Ensure correct ordering

        # Mapping dictionary to track original and new residue IDs
        mapping = {}
        count = 1
        for original_resid in original_residues:
            mapping[original_resid] = count
            pymol.cmd.alter(f"truncated_protein and resi {original_resid}", f"resi='{count}'")
            count += 1

        # Save the renumbered, truncated protein
        pymol.cmd.save(output_file, "truncated_protein")
        
        print(f"Truncated and renumbered protein saved to {output_file}")
        print("Residue mapping (original -> new):")
        for orig, new in mapping.items():
            print(f"{orig} -> {new}")
    return mapping


def truncate_pdb(pdb_file,
                 selection,
                 output_file,
                 distance_cutoff=15.0):

    with pymol2.PyMOL() as pymol:
        pymol.cmd.reinitialize()
        
        # Load the PDB file
        pymol.cmd.load(pdb_file, "protein")
        
        # Select residues within the distance cutoff
        pymol.cmd.select("close_residues", f"byres (br. ({selection}) within {distance_cutoff})")
        
        # Extract and save the selection to a temporary object
        pymol.cmd.create("truncated_protein", "close_residues")
        
        # Get unique chain IDs in the selection
        chains = []
        pymol.cmd.iterate("truncated_protein", "chains.append(chain)", space={'chains': chains})
        unique_chains = set(chains)
        
        # Mapping dictionary to track original and new residue IDs
        mapping = {}
        
        # Renumber residues in each chain separately, starting from 1 for each chain
        for chain_id in unique_chains:
            original_residues = []
            pymol.cmd.iterate(f"truncated_protein and chain {chain_id}", 
                              "original_residues.append(resi)", space={'original_residues': original_residues})
            original_residues = list(dict.fromkeys(original_residues))  # Remove duplicates, preserve order
            
            count = 1
            for original_resid in original_residues:
                mapping[(chain_id, original_resid)] = count
                pymol.cmd.alter(f"truncated_protein and chain {chain_id} and resi {original_resid}", f"resi='{count}'")
                count += 1
        
        # Save the renumbered, truncated protein
        pymol.cmd.save(output_file, "truncated_protein")
        
        print(f"Truncated and renumbered protein saved to {output_file}")
        print("Residue mapping (chain, original -> new):")
        for (chain, orig), new in mapping.items():
            print(f"{chain} {orig} -> {new}")
    return mapping


if __name__ == "__main__":

    def list_of_strings(arg):
        return arg.split(',')

    import argparse
    parser = argparse.ArgumentParser()

    parser.add_argument('-p',
                        '--pdb_file',
                        type=str,
                        help='input pdb with protein')

    parser.add_argument('-s',
                        '--selection',
                        required=True,
                        type=str,
                        help='selection in pymol language\
                            example: "chain A and resi 50" ')

    parser.add_argument('-O',
                        '--output_file',
                        type=str,
                        help='directory to output data')

    parser.add_argument('-d',
                        '--distance_cutoff',
                        type=float,
                        required=False,
                        default=15.0,
                        help='distance for cutoff')

    args = parser.parse_args()

    import os
    try:
        os.mkdir(f'{args.outdir}')
    except:
        pass

    # Usage
    residue_mapping = truncate_pdb(args.pdb_file,
                 args.selection,
                 args.output_file,
                 distance_cutoff=args.distance_cutoff)
