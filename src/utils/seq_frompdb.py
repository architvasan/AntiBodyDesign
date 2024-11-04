def get_seq_from_pdb( pdb_fn ) -> str:
    '''
    Given a pdb file, return the sequence of the protein as a string.
    '''

    to1letter = {
    "ALA":'A', "ARG":'R', "ASN":'N', "ASP":'D', "CYS":'C',
    "GLN":'Q', "GLU":'E', "GLY":'G', "HIS":'H', "ILE":'I',
    "LEU":'L', "LYS":'K', "MET":'M', "PHE":'F', "PRO":'P',
    "SER":'S', "THR":'T', "TRP":'W', "TYR":'Y', "VAL":'V' }

    seq = []
    seqstr = ''
    with open(pdb_fn) as fp:
        for line in fp:
            if line.startswith("TER"):
                seq.append(seqstr)
                seqstr = ''
            if not line.startswith("ATOM"):
                continue
            if line[12:16].strip() != "CA":
                continue
            resName = line[17:20]
            
            seqstr += to1letter[resName]

    return seq

if __name__ == "__main__":

    def list_of_strings(arg):
        return arg.split(',')

    import argparse
    parser = argparse.ArgumentParser()

    parser.add_argument('-p',
                        '--inputpdb',
                        type=str,
                        help='input pdb with protein')

    args = parser.parse_args()
    
    seq = get_seq_from_pdb(args.inputpdb)
    print(seq)
