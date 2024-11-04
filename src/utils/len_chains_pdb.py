import pymol2

def count_residues_in_chain(pdb_file,
                            chain_id):

    with pymol2.PyMOL() as pymol:
        pymol.cmd.reinitialize()
        
        # Load the PDB file
        pymol.cmd.load(pdb_file, "protein")
        
        # Create a set to store unique residue identifiers
        residues = set()
        
        # Iterate over atoms in the specified chain and add residue IDs to the set
        pymol.cmd.iterate(f"chain {chain_id}", "residues.add((chain, resv))", space={'residues': residues})
        
        # Count the number of unique residues
        num_residues = len(residues)
        
        print(f"Number of residues in chain {chain_id}: {num_residues}")
        return num_residues
