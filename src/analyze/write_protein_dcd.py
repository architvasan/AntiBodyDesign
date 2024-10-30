import MDAnalysis as mda
import numpy as np
from MDAnalysis.analysis import align
from tqdm import tqdm 
from MDAnalysis.transformations import wrap, unwrap
from MDAnalysis.transformations import center_in_box

def wrap_system(ts):
    """Manually wrap atoms to keep them inside the primary unit cell."""
    ts.positions = np.mod(ts.positions + ts.dimensions[:3], ts.dimensions[:3])
    return ts

def center_and_wrap(ts):
    # Center monomer1 in the box
    u.atoms.translate(-monomer1.center_of_mass())
    u.atoms.translate([ts.dimensions[0]/2, ts.dimensions[1]/2, ts.dimensions[2]/2])

    # Now wrap the entire system based on the center of mass
    u.atoms.positions = wrap(u.atoms.positions, ts.dimensions, compound='group')

    return ts

def extract_protein_only(pdb_inp, traj_inp, protein_only_pdb, protein_only_traj):
    u = mda.Universe(pdb_inp, traj_inp)
    protein = u.select_atoms('protein')

    with mda.Writer(protein_only_traj, protein.n_atoms) as writer:
        for ts in u.trajectory:
            writer.write(protein)

    # Optional: Write the final aligned protein structure to a PDB file
    with mda.Writer(protein_only_pdb, protein.n_atoms) as pdb_out:
        pdb_out.write(protein)
    
    print(protein_only_pdb)

    return

def align_wrap(pdb, traj, mon1_phrase, mon2_phrase, outpdb, outtraj):
    # Load the universe with your topology and trajectory files
    u = mda.Universe(pdb, traj)

    # Select only the protein
    monomer1 = u.select_atoms(mon1_phrase)
    monomer2 = u.select_atoms(mon2_phrase)

    # Create a reference structure (e.g., the first frame)
    reference = mda.Universe(pdb)
    reference_protein = reference.select_atoms('protein')
    # Select the protein
    protein = u.select_atoms("protein")
    
    # Define the transformation to center the protein in the box
    center_protein = center_in_box(protein, wrap=True)
    
    # Apply the transformation to the trajectory
    u.trajectory.add_transformations(center_protein)

    # Apply the wrapping transformation
    #wrap_transformation = wrap(u.atoms, compound='group')

    #u.trajectory.add_transformations(wrap_transformation)#wrap_system)

    # Align the trajectory to the reference structure
    aligner = align.AlignTraj(u, reference_protein, select='protein', filename=outtraj).run()

    ## Use tqdm to manually create a progress bar for the trajectory processing
    #with mda.Writer(outtraj, protein.n_atoms) as dcd_out:
    #    for ts in tqdm(u.trajectory, desc="Processing frames", total=u.trajectory.n_frames):
    #        dcd_out.write(protein)

    ## Optional: Write the final aligned protein structure to a PDB file
    # Select only the protein
    protein = u.select_atoms('protein')

    with mda.Writer(outpdb, protein.n_atoms) as pdb_out:
        pdb_out.write(protein)
    return

def main(pdb_inp, 
        traj_inp,
        mon1_phrase,
        mon2_phrase,
        outpdb,
        outtraj):

    extract_protein_only(pdb_inp, traj_inp, f'{pdb_inp}.protein.pdb', f'{traj_inp}.protein.dcd') 
    align_wrap(f'{pdb_inp}.protein.pdb', f'{traj_inp}.protein.dcd', mon1_phrase, mon2_phrase, outpdb, outtraj)
    return

if __name__ == "__main__":
    import argcomplete
    from argcomplete.completers import FilesCompleter
    import argparse
    parser = argparse.ArgumentParser(description="Process a trajectory and save only the protein.")
    parser.add_argument("-p", "--pdb_inp", help="Input PDB file")
    parser.add_argument("-t", "--traj_inp", help="Input trajectory file")
    parser.add_argument("-m1", "--m1_phrase", help="monomer 1 phrase")
    parser.add_argument("-m2", "--m2_phrase", help="monomer 2 phrase")
    parser.add_argument("-op", "--outpdb", help="Output PDB file")
    parser.add_argument("-ot", "--outtraj", help="Output trajectory file in DCD format")
    #argcomplete.autocomplete(parser)
    args = parser.parse_args()
    main(args.pdb_inp, args.traj_inp, args.m1_phrase, args.m2_phrase, args.outpdb, args.outtraj)
