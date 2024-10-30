import torch
import esm
import argparse
import os
def fold_esm(seq1,
             out_pdb,
             seq2="none",
             device=0,
            ):
    model = esm.pretrained.esmfold_v1()
    model = model.eval().cuda()
    device = torch.device(f'cuda:{device}')
    model = model.to(device)
    # Optionally, uncomment to set a chunk size for axial attention. This can help reduce memory.
    # Lower sizes will have lower memory requirements at the cost of increased speed.
    # model.set_chunk_size(128)
    
    if seq2=="none":
        sequence = seq1
    else:
        sequence = f"{seq1}:{seq2}"

    # Multimer prediction can be done with chains separated by ':'
    
    for it in range(1):
        with torch.no_grad():
            output = model.infer_pdb(sequence)
        
        with open(out_pdb, "w") as f:
            f.write(output)
        
        import biotite.structure.io as bsio
        struct = bsio.load_structure(out_pdb, extra_fields=["b_factor"])
        print(struct.b_factor.mean())  # this will be the pLDDT
        # 88.3
    return

if __name__ == "__main__": 
    parser = argparse.ArgumentParser()
    
    parser.add_argument('-T',
                        '--typeseq',
                        type=str,
                        help='how are sequences formatted? (file/str)')

    parser.add_argument('-s1',
                        '--seq1',
                        type=str,
                        help='Sequence 1')
    
    parser.add_argument('-s2',
                        '--seq2',
                        type=str,
                        required=False,
                        default="none",
                        help='Sequence 2')

    parser.add_argument('-od',
                        '--outdir',
                        type=str,
                        required=True,
                        help='output directory')

    parser.add_argument('-op',
                        '--outpdb',
                        type=str,
                        required=True,
                        help='output pdb (just file.pdb, not full path)')


    parser.add_argument('-d',
                        '--device',
                        type=str,
                        required=False,
                        default=0,
                        help='Device to place job')
    
    args = parser.parse_args()
    
    try:
        os.mkdir(f'{args.outdir}')
    except:
        print("couldn't make dir")
        pass

    if args.typeseq == 'str':
        fold_esm(args.seq1,
                 f'{args.outdir}/{args.outpdb}',
                 args.seq2,
                 args.device,
                )

    elif args.typeseq == 'file':
        with open(args.seq1, 'r') as file:
            seq1_str = file.read().replace('\n', '')

        if args.seq2 != "none":
            with open(args.seq2, 'r') as file:
                seq2_str = file.read().replace('\n', '')
        else:
            seq2_str = args.seq2

        fold_esm(seq1_str,
         f'{args.outdir}/{args.outpdb}',
         seq2_str,
         args.device,
        )

