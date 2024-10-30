import pandas as pd

def id_cdrloop(data,
                rowit,
                cdrloop,
                chaintarget):
    
    data_it = data.loc[rowit]
    data_it_chain = data_it[chaintarget]
    data_it_patt = data_it[cdrloop]
    rid_init = data_it_chain.find(data_it_patt)
    rid_fin = rid_init + len(data_it_patt) - 1
    return rid_init, rid_fin

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-i',
                        '--inputfil',
                        type=str,
                        help='input file')

    parser.add_argument('-r',
                        '--rowit',
                        type=int,
                        help='row index')

    parser.add_argument('-c',
                        '--cdrloop',
                        type=str,
                        help='column for cdrloop\
                              options: \
                                "heavy_cdr1"\
                                "heavy_cdr2"\
                                "heavy_cdr3"\
                                "light_cdr1"\
                                "light_cdr2"\
                                "light_cdr3" ')

    parser.add_argument('-T',
                        '--chaintarget',
                        type=str,
                        help='column for chaintarget \
                              options: \
                                "heavy_chain" \
                                "light_chain"\
                                ')


    args = parser.parse_args()
    data = pd.read_csv(args.inputfil)

    rid_init, rid_fin = id_cdrloop(data, args.rowit, args.cdrloop, args.chaintarget)
    print(f"Resids: {rid_init}:{rid_fin}")
