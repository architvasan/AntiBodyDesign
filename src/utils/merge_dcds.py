import MDAnalysis as mda
from MDAnalysis.coordinates.DCD import DCDWriter

def main(inputpdb, inputtraj_list, outtraj):
    # Load the first universe to get the topology
    u = mda.Universe(inputpdb, inputtraj_list)
    
    # Open a new DCDWriter for the output file
    with DCDWriter(outtraj, u.atoms.n_atoms) as writer:
        # Iterate over each frame in the current DCD file
        for ts in u.trajectory:
            # Write the current frame to the new DCD file
            writer.write(u.atoms)
    


if __name__ == "__main__":

    def list_of_strings(arg):
        return arg.split(',')

    import argparse
    parser = argparse.ArgumentParser()

    parser.add_argument('-p',
                        '--inputpdb',
                        type=str,
                        help='input pdb with protein')

    parser.add_argument('-T',
                        '--inputtraj_list',
                        required=False,
                        type=list_of_strings,
                        default=['none', 'none'],
                        help='trajectory file list (dcd format)')

    parser.add_argument('-O',
                        '--outdir',
                        type=str,
                        help='directory to output data')

    parser.add_argument('-o',
                        '--out',
                        type=str,
                        help='file to store merged dcd')

    args = parser.parse_args()

    import os
    try:
        os.mkdir(f'{args.outdir}')
    except:
        pass

    main(args.inputpdb,
         args.inputtraj_list,
         f'{args.outdir}/{args.out}',
        )
