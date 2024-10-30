import glob

import pandas as pd
import residue_analysis as ra

def main(inpdir):
    pdb_list = glob.glob(f'{inpdir}/peptides*dldesign*pdb')
    res_analysis = {'Positives': [],
                    'Negatives': [],
                    'Polars': [],
                    'NonPolars': []}
    for pdbfil in pdb_list:
        pos, neg, pol, nonpol = ra.residuecount(pdbfil)
        res_analysis['Positives'].append(pos)
        res_analysis['Negatives'].append(neg)
        res_analysis['Polars'].append(pol)
        res_analysis['NonPolars'].append(nonpol)

    return res_analysis

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
   
    res_analysis = main(args.inpdir)

    df_res_analysis = pd.DataFrame(res_analysis)

    df_res_analysis.to_csv(f'{args.outdir}/{args.outfile}')



