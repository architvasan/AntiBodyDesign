import MDAnalysis as mda
from MDAnalysis.coordinates.memory import MemoryReader
import numpy as np
def main(inppdb,
        inptraj,
        framefile,
        outdir,
        outtraj):

    # Load the original trajectory
    print(inppdb)
    print(inptraj)
    u = mda.Universe(inppdb, inptraj)
    
    frame_indices = [int(i) for i in np.loadtxt(framefile)]
    
    # Extract and store selected frames
    selected_frames = []
    for frame in frame_indices:
        u.trajectory[frame]
        selected_frames.append(u.atoms.positions.copy())
    
    # Convert list of frames to a single numpy array with shape (n_frames, n_atoms, 3)
    combined_positions = np.array(selected_frames)
    
    # Create a new trajectory file with the selected frames
    with mda.Writer(f'{outdir}/{outtraj}', n_atoms=u.atoms.n_atoms) as W:
        for positions in combined_positions:
            u.atoms.positions = positions  # Update the positions
            W.write(u.atoms)


if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser()

    parser.add_argument('-p',
                        '--pdbinp',
                        type=str,
                        help='input pdb')

    parser.add_argument('-t',
                        '--trajinp',
                        type=str,
                        help='input trajectory (dcd)')

    parser.add_argument('-F',
                        '--framefile',
                        type=str,
                        help='file with list of frames')

    parser.add_argument('-O',
                        '--outdir',
                        type=str,
                        help='directory for output trajs')

    parser.add_argument('-ot',
                        '--outtraj',
                        type=str,
                        help='output trajectory')

    args = parser.parse_args()

    try:
        import os
        os.mkdir(args.outdir)
    except:
        pass

    main(args.pdbinp,
         args.trajinp,
         args.framefile,
         args.outdir,
         args.outtraj)
