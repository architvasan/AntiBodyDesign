import MDAnalysis as mda
from MDAnalysis.lib.mdamath import make_whole
from MDAnalysis.lib.util import unique_int_1d
from MDAnalysis.analysis.distances import distance_array
import numpy as np
from tqdm import tqdm

def calc_cont(u, selection1, selection2, cutoff=3.5):

    # Select the two chains by their chain identifiers (e.g., 'A' and 'B')
    chain_A = u.select_atoms(selection1)
    chain_B = u.select_atoms(selection2)
    res_A = sorted(list(set([res.resid for res in chain_A])))
    res_B = sorted(list(set([res.resid for res in chain_B])))
    #res_conts_A = np.zeros(len(res_A)+1)#{}
    #res_conts_B = np.zeros(len(res_B)+1)#{}
    #cont_A_B_res = np.zeros((len(res_A)+1, len(res_B)+1))

    conts_per_frame = []
    A_start = res_A[0]
    B_start = res_B[0]
    del chain_B

    for ts in (u.trajectory):
        contacts_total_ts = 0
        chain_B = u.select_atoms(f"{selection2} and around {cutoff} ({selection1})") 
        # Calculate the pairwise distances between all atoms in chain A and chain B
        distances = distance_array(chain_A.positions, chain_B.positions)
        # Find the pairs of atoms with a distance below a certain threshold (e.g., 4.0 Å)
        contacts = (distances < cutoff)
        conts_per_frame.append(np.sum(contacts)) 
        continue
        #res_conts_A_ts = np.zeros(len(res_A)+1)
        cont_A_B_res_ts = np.zeros((len(res_A)+1, len(res_B)+1))
        # Print or analyze the contacts
        for i, atom_A in enumerate(chain_A):
            for j, atom_B in enumerate(chain_B):
                if contacts[i, j]:
                    #print(f"Contact between {atom_A.resname}{atom_A.resid} (chain A) and {atom_B.resname}{atom_B.resid} (chain B) at distance {distances[i, j]:.2f} Å")
                    A_rid_it = int(atom_A.resid)-A_start
                    B_rid_it = int(atom_B.resid)-B_start

                    if cont_A_B_res_ts[A_rid_it][B_rid_it]==0:
                        contacts_total_ts+=1
        del chain_B
        conts_per_frame.append(contacts_total_ts)
    
    print(conts_per_frame)
    return conts_per_frame

def main(inputpdb, 
         inputtraj_list,
         sel_A,
         sel_B,
         threshold,
         out_contacts_per_frame,
         ):

    # Load the PDB file
    u = mda.Universe(inputpdb, inputtraj_list)
    
    conts_per_frame = calc_cont(u, sel_A, sel_B, threshold)
    # Select the two chains by their chain identifiers (e.g., 'A' and 'B')
    
    #np.savetxt(out_contacts_per_frame, conts_per_frame)


    return conts_per_frame

if __name__ == "__main__":

    def list_of_strings(arg):
        return arg.split(',')

    import argparse
    parser = argparse.ArgumentParser()

    parser.add_argument('-p',
                        '--inputpdb',
                        type=str,
                        help='input pdb with protein')

#    parser.add_argument('-t',
#                        '--inputtraj',
#                        required=False,
#                        type=str,
#                        default='none',
#                        help='single trajectory file')

    parser.add_argument('-T',
                        '--inputtraj_list',
                        required=False,
                        type=list_of_strings,
                        default=['none', 'none'],
                        help='trajectory file list (dcd format)')

    parser.add_argument('-sA',
                        '--selA',
                        type=str,
                        help='phrase for seelction A (in mdanalysis language)')

    parser.add_argument('-sB',
                        '--selB',
                        type=str,
                        help='phrase for seelction A (in mdanalysis language)')

    parser.add_argument('-c',
                        '--cutoff',
                        type=float,
                        help='cutoff for judging a contact or not (3.5 for heavy atoms)')

    parser.add_argument('-O',
                        '--outdir',
                        type=str,
                        help='directory to output data')

    parser.add_argument('-of',
                        '--outframes',
                        type=str,
                        help='file to output contacts per frame (dont use path, just file name)')

    args = parser.parse_args()

    import os
    try:
        os.mkdir(f'{args.outdir}')
    except:
        pass

    main(args.inputpdb,
         args.inputtraj_list,
         args.selA,
         args.selB,
         args.cutoff,
         f'{args.outdir}/{args.outframes}',
         )

