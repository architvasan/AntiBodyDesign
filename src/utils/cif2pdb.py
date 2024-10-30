import pymol2
import os

def remove_extension(filename):
    return os.path.splitext(filename)[0]

def cif2pdb(filename): 
    with pymol2.PyMOL() as pymol:
         pymol.cmd.load(filename, 'myprotein')
         pymol.cmd.save(filename.replace('.cif', '.pdb'), selection='myprotein')

if __name__ =="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-i',
                        '--inputfil',
                        type=str,
                        help='input cif')
    args = parser.parse_args()
    cif2pdb(args.inputfil)
