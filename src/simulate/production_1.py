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

def prod1(inpcrd_fil,
            prmtop_fil,
            d_ind,
            eq_chkpt,
            prod_steps,
            prod_dcd,
            prod_rst,
            prod_chkpt,
            prod_st,
            prod_log):

    system, prmtop, inpcrd = load_amber_files(inpcrd_fil,
                                            prmtop_fil)
    
    eq_simulation, integrator = setup_sim_nomin(system,
                                                prmtop,
                                                inpcrd,
                                                d_ind=d_ind)
    
    barostat = system.addForce(MonteCarloBarostat(1*atmosphere, 300*kelvin))
    
    prod_1_sim = run_prod(eq_simulation, 
                            prod_steps,
                            eq_chkpt,
                            prod_dcd,
                            prod_rst,
                            prod_chkpt,
                            prod_st,
                            prod_log)   
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
                                                 
    parser.add_argument('-n',
                        '--numsteps',
                        type=int,
                        help='Number of md steps to run')
    
    args = parser.parse_args()
    
    try:
        os.mkdir(f'{args.rep_loc}/prod1')
    except:
        pass
    inpcrd_fil = f'{args.struct}/inpcrd1'
    prmtop_fil = f'{args.struct}/prmtop1'
    d_ind = args.device
    eq_chkpt = f'{args.rep_loc}/prod0/prod0.rst.chk'
    prod_steps = args.numsteps #100000000
    prod_dcd = f'{args.rep_loc}/prod1/prod1.dcd'
    prod_rst = f'{args.rep_loc}/prod1/prod1.rst.chk'
    prod_chkpt = f'{args.rep_loc}/prod1/prod1.chk'
    prod_st = f'{args.rep_loc}/prod1/prod1.state'
    prod_log = f'{args.rep_loc}/prod1/prod1.csv'

    prod1(inpcrd_fil,
      prmtop_fil,
      d_ind,
      eq_chkpt,
      prod_steps,
      prod_dcd,
      prod_rst,
      prod_chkpt,
      prod_st,
      prod_log)




























if False:
    try:
        os.mkdir('prod1')
    except:
        pass
    inpcrd_fil = '../struct/inpcrd1'
    prmtop_fil = '../struct/prmtop1'
    #eq_st = 'prod0/prod0.state'
    eq_chkpt = 'prod0/prod0.rst.chk'
    d_ind = '2'
    prod_steps = 100000000
    prod_chkpt = 'prod1/prod1.chk'
    prod_st = 'prod1/prod1.state'
    prod_dcd = 'prod1/prod1.dcd'
    prod_rst = 'prod1/prod1.rst.chk'
    prod_log = 'prod1/prod1.csv'
    
    system, prmtop, inpcrd = load_amber_files(inpcrd_fil,
                                            prmtop_fil)
    
    eq_simulation, integrator = setup_sim_nomin(system,
                                                prmtop,
                                                inpcrd,
                                                d_ind=d_ind)
    
    eq_simulation.context.setParameter('k', 0)
    
    prod_1_sim = run_prod(eq_simulation, 
                            prod_steps,
                            eq_chkpt,
                            prod_dcd,
                            prod_rst,
                            prod_chkpt,
                            prod_st,
                            prod_log)
