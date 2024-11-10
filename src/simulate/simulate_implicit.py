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
from simulation_funcs import *
import argparse

def running(input_pdb, simulation_time, output_dcd, d_ind): 
    system, pdb, forcefield = system_implicit(input_pdb)
    simulation, potential_energy = sim_implicit(
                                        system,
                                        pdb,
                                        simulation_time,
                                        output_dcd,
                                        d_ind,
                                        )
    return


if __name__ == "__main__": 
    parser = argparse.ArgumentParser()
    
    parser.add_argument('-p',
                        '--pdb_inp',
                        type=str,
                        help='Location of input pdb')
    
    parser.add_argument('-t',
                        '--time_sim',
                        type=int,
                        required=False,
                        default=500000,
                        help='Simulation time (500,000: 1ns)')

    parser.add_argument('-od',
                        '--outputdcd',
                        type=str,
                        help='outputdcd')

    parser.add_argument('-d',
                        '--device',
                        type=str,
                        help='Device to place job')
    
    args = parser.parse_args()
    
    running(args.pdb_inp, 
            args.time_sim, 
            args.outputdcd,
            args.device)