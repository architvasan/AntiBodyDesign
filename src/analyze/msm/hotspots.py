import pandas as pd

def main(datainp, cutoff, resid_col, prob_col, output):
    df = pd.read_csv(datainp)
    df_top = df[df[prob_col]>cutoff]
    df_top.to_csv(output, index=False)
    return df

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()

    parser.add_argument('-d',
                        '--datainp',
                        type=str,
                        help='datafile from calculating residue contact probabilities')

    parser.add_argument('-c',
                        '--cutoff',
                        type=float,
                        help='prob cutoff')

    parser.add_argument('-r',
                        '--rcol',
                        type=str,
                        help='resid column')

    parser.add_argument('-p',
                        '--pcol',
                        type=str,
                        help='prob column')

    parser.add_argument('-O',
                        '--outdir',
                        type=str,
                        help='directory for output data')

    parser.add_argument('-o',
                        '--outfile',
                        type=str,
                        help='output file')

    args = parser.parse_args()

    try:
        os.mkdir(f'{args.outdir}')
    except:
        pass

    try:
        os.mkdir(f'{args.imdir}')
    except:
        pass

    main(args.datainp,
         args.cutoff,
         args.rcol,
         args.pcol,
         f'{args.outdir}/{args.outfile}',
         )
