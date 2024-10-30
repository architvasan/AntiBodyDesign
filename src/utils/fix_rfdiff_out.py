import MDAnalysis as mda
import os

import pymol2
from pymol import cmd

def modify_chain(input_file, resid_cutoff):
    with pymol2.PyMOL() as pymol:
        # Load the input PDB file
        print(input_file)
        pymol.cmd.load(input_file, "structure")
        total_residues = pymol.cmd.count_atoms("(name CA)")
        # Select the target residues in chain A (residues 227-468)
        pymol.cmd.select(f"target_residues", f"chain A and resi {resid_cutoff}-{total_residues}")
        print(f"chain A and resi {resid_cutoff}-{total_residues}")

        # Change the chain ID to B for the selected residues
        pymol.cmd.alter("target_residues", "chain = 'B'")

        # Rebuild to apply the changes
        pymol.cmd.rebuild()

        # Save the modified structure to the output file
        pymol.cmd.save(input_file, "structure")


