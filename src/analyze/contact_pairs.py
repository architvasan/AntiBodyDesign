import MDAnalysis as mda
from MDAnalysis.lib.mdamath import make_whole
from MDAnalysis.lib.util import unique_int_1d
from MDAnalysis.analysis.distances import distance_array
import numpy as np
from tqdm import tqdm
import pandas as pd
from collections import Counter

def pairs_to_frequency_df(pairs):
    # Count the frequency of each pair
    pair_counts = Counter(pairs)
    
    # Create a DataFrame
    df = pd.DataFrame(pair_counts.items(), columns=["Pair", "Frequency"])
    
    # Optionally sort the DataFrame by frequency (descending)
    df = df.sort_values(by="Frequency", ascending=False).reset_index(drop=True)
    
    return df

def three_to_one(three_letter_code):
    # Mapping of three-letter to one-letter amino acid codes
    aa_dict = {
        "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
        "GLU": "E", "GLN": "Q", "GLY": "G", "HIS": "H", "ILE": "I",
        "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
        "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
        # Special cases
        "HIE": "H", "HID": "H", "SEC": "U", "PYL": "O", "ASX": "B",
        "GLX": "Z", "XLE": "J", "TER": "*",
    }
    
    # Convert the input code to uppercase and return the one-letter code
    return aa_dict.get(three_letter_code.upper(), "Unknown")

def calc_cont_pairs(inputpdb,
                    selection1,
                    selection2, 
                    chainA, 
                    chainB, 
                    cutoff=3.5):

    dict_pairs = {'chainA': [],
                  'res_idxA': [],
                  'chainB': [],
                  'res_idxB': [],
                  'connection_type': [],
                  'confidence': [],
                  'min_distance_angstrom': [],
                  'max_distance_angstrom': [],
                  'comment': [],
                  'restraint_id': [],
                }

    u = mda.Universe(inputpdb)
    # Select the two chains by their chain identifiers (e.g., 'A' and 'B')
    chain_A = u.select_atoms(selection1)
    chain_B_total = u.select_atoms(selection2)
    res_A = sorted(list(set([res.resid for res in chain_A])))
    res_B = sorted(list(set([res.resid for res in chain_B_total])))
    A_start = res_A[0]
    B_start = res_B[0]
    del chain_B_total
    contact_pairs = []

    chain_B = u.select_atoms(f"{selection2} and around {cutoff} ({selection1})") 
    # Calculate the pairwise distances between all atoms in chain A and chain B
    distances = distance_array(chain_A.positions, chain_B.positions)
    # Find the pairs of atoms with a distance below a certain threshold (e.g., 4.0 Å)
    contacts = (distances < cutoff)
    # Print or analyze the contacts
    rest_id = 0
    for i, atom_A in enumerate(chain_A):
        for j, atom_B in enumerate(chain_B):
            if contacts[i, j]:
                A_rid_it = int(atom_A.resid)
                B_rid_it = int(atom_B.resid)
                A_rnm_trip = atom_A.residue.resname
                B_rnm_trip = atom_B.residue.resname
                
                A_rnm_sing = three_to_one(A_rnm_trip)
                B_rnm_sing = three_to_one(B_rnm_trip)
                A_atom_name = atom_A.name
                B_atom_name = atom_B.name
                if abs(A_rid_it - B_rid_it)>4:
                    A_sel_it = u.select_atoms(f'resid {A_rid_it} and name CA')
                    B_sel_it = u.select_atoms(f'resid {B_rid_it} and name CA')
                    dist = distance_array(A_sel_it.positions, B_sel_it.positions)[0][0] 
                    print(dist)
                    del(A_sel_it)
                    del(B_sel_it)
                    pair_it = (f'{A_rnm_sing}{A_rid_it}',
                               f'{B_rnm_sing}{B_rid_it}')
                    if pair_it not in contact_pairs:
                        rest_id+=1
                        contact_pairs.append((f'{A_rnm_sing}{A_rid_it}',
                                           f'{B_rnm_sing}{B_rid_it}'))
                        dict_pairs["chainA"].append(chainA)
                        dict_pairs["res_idxA"].append(f'{A_rnm_sing}{A_rid_it-1}')
                        dict_pairs["chainB"].append(chainB)
                        dict_pairs["res_idxB"].append(f'{B_rnm_sing}{B_rid_it-1}')
                        dict_pairs["connection_type"].append('contact')
                        dict_pairs["confidence"].append(1.0)
                        dict_pairs["min_distance_angstrom"].append(0.0)
                        dict_pairs["max_distance_angstrom"].append(dist)
                        dict_pairs["comment"].append(f'{A_rnm_sing}{A_rid_it}:{B_rnm_sing}{B_rid_it}'\
                                                     f':{cutoff}')
                        dict_pairs["restraint_id"].append(f'restraint{rest_id}')
                    else:
                        contact_pairs.append((f'{A_rnm_sing}{A_rid_it}',
                                           f'{B_rnm_sing}{B_rid_it}'))
    del chain_B
    freq_pairs_df = pairs_to_frequency_df(contact_pairs)
    df_pairs_rests = pd.DataFrame(dict_pairs)

    return freq_pairs_df, df_pairs_rests

if __name__ == "__main__":

    def list_of_strings(arg):
        return arg.split(',')

    import argparse
    parser = argparse.ArgumentParser()

    parser.add_argument('-p',
                        '--inputpdb',
                        type=str,
                        help='input pdb with protein')
    
    parser.add_argument('-sA',
                        '--selA',
                        type=str,
                        help='phrase for seelction A (in mdanalysis language)')

    parser.add_argument('-sB',
                        '--selB',
                        type=str,
                        help='phrase for seelction B (in mdanalysis language)')

    parser.add_argument('-CA',
                        '--chainA',
                        type=str,
                        help='chain for seelction A')

    parser.add_argument('-CB',
                        '--chainB',
                        type=str,
                        help='chain for seelction B')

    parser.add_argument('-c',
                     '--cutoff',
                     type=float,
                     required=False,
                     default=5.5,
                     help='cutoff (default: 5.5)')

    parser.add_argument('-of',
                        '--outfreq',
                        type=str,
                        help='pair freq data outfile')

    parser.add_argument('-or',
                        '--outrest',
                        type=str,
                        help='restraints format outfile?')

    args = parser.parse_args()

    freq_pairs_df,df_pairs_rests = calc_cont_pairs(args.inputpdb,
                                                   args.selA, 
                                                   args.selB, 
                                                   chainA=args.chainA,
                                                   chainB=args.chainB,
                                                   cutoff=args.cutoff)
    freq_pairs_df.to_csv(args.outfreq, index=False)
    df_pairs_rests.to_csv(args.outrest, index=False)
    # Write the pairs to the file
    # with open(args.outfile, "w") as file:
        # for pair in int_pairs:
            # Convert each pair to a string and write it to the file
            # file.write(f"{pair[0]}, {pair[1]}\n")