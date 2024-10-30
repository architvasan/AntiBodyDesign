import MDAnalysis as mda
import time
from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout, exit, stderr
from parmed import unit as u
from copy import deepcopy
import sys
from sys import stdout
#from openff.toolkit import Molecule
#from openmmforcefields.generators import GAFFTemplateGenerator
import pandas as pd
import numpy as np
from parmed import load_file, unit as u
from .simulation_funcs import *
import argparse

def running_eq(inpcrd_fil, prmtop_fil, eq_st, eq_chkpt, d_ind): 
    eq_simulation = run_eq(inpcrd_fil, prmtop_fil, eq_st, eq_chkpt, d_ind=d_ind)
    return


if __name__ == "__main__": 
    parser = argparse.ArgumentParser()
    
    parser.add_argument('-R',
                        '--rep_loc',
                        type=str,
                        help='Directory for replica')
    
    parser.add_argument('-s',
                        '--struct',
                        type=str,
                        help='Directory for structural files')


    parser.add_argument('-d',
                        '--device',
                        type=str,
                        help='Device to place job')
    
    args = parser.parse_args()
    
    try:
        os.mkdir(f'{args.rep_loc}/prod0')
    except:
        pass

    inpcrd_fil = f'{args.struc}/inpcrd1'
    prmtop_fil = f'{args.struct}/prmtop1'
    try:
        os.mkdir(f'{args.rep_loc}/eq')
    except:
        pass

    eq_st = f'{args.rep_loc}/eq/eq.state'
    eq_chkpt = f'{args.rep_loc}/eq/eq.chk'
    eq_simulation = run_eq(inpcrd_fil, prmtop_fil, eq_st, eq_chkpt, d_ind= args.device)
