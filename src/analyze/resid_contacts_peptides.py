import MDAnalysis as mda
from MDAnalysis.lib.mdamath import make_whole
from MDAnalysis.lib.util import unique_int_1d
from MDAnalysis.analysis.distances import distance_array
import glob
import pandas as pd
import contacts_A_B as contsab
from tqdm import tqdm

def main(inpdir):
    pdb_list = glob.glob(f'{inpdir}/peptides*dldesign*pdb')
    u = mda.Universe(pdb_list[0], pdb_list[0])
    A_rids_grnd = [it for it in range(1, 308)]
    A_rid_probs, B_rid_probs, cont_mat = contsab.calc_rid_cont(u, 'segid B', 'segid A')
    print(list(A_rid_probs.values()))
    A_rid_probs_tot = {rid: prob for rid, prob in zip(A_rids_grnd, list(A_rid_probs.values()))}

    for pdbfil in tqdm(pdb_list):
        u = mda.Universe(pdbfil, pdbfil)
        A_rid_probs, B_rid_probs_it, cont_mat_it = contsab.calc_rid_cont(u, 'segid B', 'segid A')

        A_rid_probs_fixed = {rid: prob for rid, prob in zip(A_rids_grnd, list(A_rid_probs.values()))}

        A_rid_probs = {key: A_rid_probs_fixed.get(key, 0) + A_rid_probs_fixed.get(key, 0) for key in set(A_rid_probs_fixed) | set(A_rid_probs_fixed)}

    A_rid_probs_normed = {key: value / len(pdb_list) for key, value in A_rid_probs.items()}

    return A_rid_probs_normed

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()

    parser.add_argument('-p',
                        '--inpdir',
                        type=str,
                        help='input directory with peptides')

    parser.add_argument('-O',
                        '--outdir',
                        type=str,
                        help='directory to output data')

    parser.add_argument('-o',
                        '--outfile',
                        type=str,
                        help='file to output')


    args = parser.parse_args()

    A_rid_probs = main(args.inpdir)

    with open(f'{args.outdir}/{args.outfile}', 'w') as file:
        file.write('resid,contactprob\n')
        for resid, prob in A_rid_probs.items():
            file.write(f'{resid},{prob:.4f}\n')
