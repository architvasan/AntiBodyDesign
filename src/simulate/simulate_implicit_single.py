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
import simulation_funcs as sfuncs
#from simulation_funcs import *
import argparse
from tqdm import tqdm
import os

def running(input_pdb, simulation_time, output_dcd, d_ind): 
    system, pdb, forcefield = sfuncs.system_implicit(input_pdb)
    simulation, potential_energy = sfuncs.sim_implicit(
                                        system,
                                        pdb,
                                        simulation_time,
                                        output_dcd,
                                        f'{output_dcd}.pdb',
                                        f'{output_dcd}.log',
                                        d_ind
                                        )
    return

def simulate_struct(
                    input_pdb,
                    simulation_time,
                    output_dcd,
                    output_pdb,
                    output_log,
                    d_ind,
                    ):

    system, pdb, forcefield = sfuncs.system_implicit(input_pdb)
    simulation, potential_energy = sfuncs.sim_implicit(
                                        system,
                                        pdb,
                                        simulation_time,
                                        output_dcd,
                                        output_pdb,
                                        output_log,
                                        d_ind,
                                        )

    return potential_energy

def simulate_pipeline(
                     df_paths,
                     output_dir,
                     simulation_time=250000,
                     gpu_per_node=4,
                     use_mpi=False,
                     device="0",
                      ):
    '''
    1. load in df of paths for all pdbs to simulate
    2. Split list of "fixed" pdbs across multi ranks
    3. Iterate over list of pdbs on each rank
      and simulate each structure
    4. Write to logfile potential energy
    5. Figure out how to extract cdr loop info from each structure. 
        Can just go to directory, get cdr loop from .csv file
        Then find the rid and rfin for that cdr in the original chai file
        And use this rid and rfin to identify the cdr in the new file
    '''
    if use_mpi:
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        size = comm.Get_size()
        rank = comm.Get_rank()
        device = str(rank % gpu_per_node)

        df_paths_chunks = np.array_split(df_paths,
                                         size)
        df_paths_rank = df_paths_chunks[rank]
        print(df_paths_rank)
    else:
        #os.environ['CUDA_VISIBLE_DEVICES']=str(device)
        df_paths_rank = df_paths
        print(device)

    structures = []
    potential_energies = []
    for rowit in tqdm(range(len(df_paths_rank))):
        df_it = df_paths_rank.iloc[rowit]
        output_name_it = df_it['VariableName']
        pdb_path_it = df_it['Location']
        out_dcd = f'{output_dir}/{output_name_it}.dcd'
        out_pdb = f'{output_dir}/{output_name_it}.pdb'
        out_log = f'{output_dir}/{output_name_it}.log'
        pot_en = simulate_struct(
                    pdb_path_it,
                    simulation_time,
                    out_dcd,
                    out_pdb,
                    out_log,
                    device,
                    )
        structures.append(out_pdb)
        potential_energies.append(pot_en)
    
    return structures, potential_energies



if True:
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
