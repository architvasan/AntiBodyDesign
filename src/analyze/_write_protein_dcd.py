import MDAnalysis as mda
import argparse
from MDAnalysis.transformations import wrap#, align_along_axis
from MDAnalysis.analysis import align
from tqdm import tqdm

# Define the wrap transformation
def wrap_system(ts):
    ts.positions = wrap(ts.positions, ts.dimensions)
    return ts

def main(pdb_inp, traj_inp, outpdb, outtraj):
    # Load the universe with your topology and trajectory files
    u = mda.Universe(pdb_inp, traj_inp)
    
    # Select only the protein
    protein = u.select_atoms('protein')
    reference = u.select_atoms('protein').copy()


    # Wrap the entire system to the primary unit cell
    #wrap_transform = wrap(u.trajectory.ts.dimensions)
    
    # Apply the transformations
    u.trajectory.add_transformations(wrap_system)

    # Perform the alignment and save the trajectory
    aligner = align.AlignTraj(u, reference, select='protein', in_memory=True).run()
    
    # Write the protein-only PDB file
    with mda.Writer(outpdb, protein.n_atoms) as pdb_out:
        pdb_out.write(protein)
    
    # Write the protein-only trajectory as a DCD file
    with mda.Writer(outtraj, protein.n_atoms) as dcd_out:
        for ts in tqdm(u.trajectory):
            dcd_out.write(protein)
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    
    parser.add_argument('-p',
                        '--pdb_inp',
                        type=str,
                        help='input pdb with whole system')
    
    parser.add_argument('-t',
                        '--traj_inp',
                        type=str,
                        help='trajectory file (dcd format)')
                                                               
                                                               
    parser.add_argument('-op',
                        '--outpdb',
                        type=str,
                        help='output pdb')
                                                 
    parser.add_argument('-ot',
                        '--outtraj',
                        type=str,
                        help='output traj (dcd)')
    
    args = parser.parse_args()
    

    main(args.pdb_inp, args.traj_inp, args.outpdb, args.outtraj)
        
