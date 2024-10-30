import MDAnalysis as mda
from MDAnalysis.lib.mdamath import make_whole
from MDAnalysis.lib.util import unique_int_1d
from MDAnalysis.analysis.distances import distance_array
import numpy as np
from tqdm import tqdm

def calc_rid_cont(u, selection1, selection2, cutoff=3.5):

    # Select the two chains by their chain identifiers (e.g., 'A' and 'B')
    chain_A = u.select_atoms(selection1)
    chain_B = u.select_atoms(selection2)
    res_A = sorted(list(set([res.resid for res in chain_A])))
    print(res_A)
    res_B = sorted(list(set([res.resid for res in chain_B])))
    print(res_B)
    res_conts_A = np.zeros(len(res_A)+1)#{}
    res_conts_B = np.zeros(len(res_B)+1)#{}
    cont_A_B_res = np.zeros((len(res_A)+1, len(res_B)+1))

    A_start = res_A[0]
    B_start = res_B[0]
    del chain_B

    for ts in tqdm(u.trajectory[::10]):
        chain_B = u.select_atoms(f"{selection2} and around {cutoff} ({selection1})") 
        # Calculate the pairwise distances between all atoms in chain A and chain B
        distances = distance_array(chain_A.positions, chain_B.positions)
        # Find the pairs of atoms with a distance below a certain threshold (e.g., 4.0 Å)
        contacts = (distances < cutoff)
    
        cont_A_B_res_ts = np.zeros((len(res_A)+1, len(res_B)+1))
        # Print or analyze the contacts
        for i, atom_A in enumerate(chain_A):
            for j, atom_B in enumerate(chain_B):
                if contacts[i, j]:
                    #print(f"Contact between {atom_A.resname}{atom_A.resid} (chain A) and {atom_B.resname}{atom_B.resid} (chain B) at distance {distances[i, j]:.2f} Å")
                    A_rid_it = int(atom_A.resid)-A_start
                    B_rid_it = int(atom_B.resid)-B_start

                    if cont_A_B_res_ts[A_rid_it][B_rid_it]==0:
                        res_conts_A[A_rid_it] += 1
                        res_conts_B[B_rid_it] += 1
                        cont_A_B_res[A_rid_it][B_rid_it]+=1
                        cont_A_B_res_ts[A_rid_it][B_rid_it]+=1
        del chain_B
    
    res_conts_A /= len(u.trajectory[::10])
    res_conts_B /= len(u.trajectory[::10])
    cont_A_B_res /= len(u.trajectory[::10])

    # Calculate contact probabilities
    A_contact_probs = {resid: prob for resid, prob in zip(res_A, res_conts_A)}
    B_contact_probs = {resid: prob for resid, prob in zip(res_B, res_conts_B)}

    return A_contact_probs, B_contact_probs, cont_A_B_res



def _calculate_residue_contacts(u, selection1, selection2, cutoff=4.5):
    # Define two groups of residues based on selections
    group1 = u.select_atoms(selection1).residues
    group2 = u.select_atoms(selection2).residues
    
    print(group1)
    print(group2)
    # Initialize an empty contact map
    contact_map = np.zeros((len(group1), len(group2)))
    del(group2)
    print(np.shape(contact_map))
    # Loop over all frames in the trajectory
    for ts in tqdm(u.trajectory):
        group2 = u.select_atoms(f"{selection2} and around {cutoff} ({selection1})") 
        if len(group2)==0:
            continue
        for i, res1 in enumerate(group1):
            for j, res2 in enumerate(group2):
                # Check if any heavy atoms between the two residues are within the cutoff
                if np.any(np.linalg.norm(res1.positions[:, None] - res2.positions, axis=-1) < cutoff):
                    contact_map[int(res1.resid), int(res2.resid)] += 1
        del group2
    # Normalize the contact map by the number of frames
    contact_map /= len(u.trajectory)
    
    print(contact_map)
    return contact_map, group1, group2

def main(inputpdb, 
         inputtraj_list,
         sel_A,
         sel_B,
         threshold,
         out_contact_map,
         out_A_probs,
         out_B_probs):

    # Load the PDB file
    u = mda.Universe(inputpdb, inputtraj_list)
    
    A_rid_probs, B_rid_probs, cont_A_B_res = calc_rid_cont(u, sel_A, sel_B, threshold)
    # Select the two chains by their chain identifiers (e.g., 'A' and 'B')
    
    np.savetxt(out_contact_map, cont_A_B_res)

    # Write results to a file
    with open(out_A_probs, 'w') as file:
        file.write('resid #, contact prob\n')
        for resid, prob in A_rid_probs.items():
            file.write(f'{resid}, {prob:.4f}\n')

        # Write results to a file
    with open(out_B_probs, 'w') as file:
        file.write('resid #, contact prob\n')
        for resid, prob in B_rid_probs.items():
            file.write(f'{resid}, {prob:.4f}\n')

    return cont_A_B_res, A_rid_probs, B_rid_probs


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

    parser.add_argument('-oc',
                        '--outc_map',
                        type=str,
                        help='file to output contactmap (dont use path, just file name)')

    parser.add_argument('-oa',
                        '--outa_prob',
                        type=str,
                        help='file to output selection A resid probs (dont use path, just file name)')

    parser.add_argument('-ob',
                        '--outb_prob',
                        type=str,
                        help='file to output selection B resid probs (dont use path, just file name)')

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
         f'{args.outdir}/{args.outc_map}',
         f'{args.outdir}/{args.outa_prob}',
         f'{args.outdir}/{args.outb_prob}')

