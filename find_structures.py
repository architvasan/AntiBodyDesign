import os
import pandas as pd
import json
import numpy as np
from tqdm import tqdm
#import src.simulate.simulation_funcs as sfuncs
import os
#from mpi4py import MPI

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



def find_fixed_pdb_and_generate_variable(pattern,
                                         outpatt,
                                         directory='trials/T14_ant3HFM_body4NCO'):
    result = {}
    
    # Traverse the directory structure
    for root, dirs, files in os.walk(directory):
        # Check if 'fixed.pdb' exists in the current directory
        if pattern in files:
            # Get the relative path from the base directory
            relative_path = os.path.relpath(root, directory)
            
            # Convert the relative path to a variable format
            path_parts = relative_path.split(os.sep)
            variable_name = f"{outpatt}_{'_'.join(path_parts)}_output"
            
            # Save the variable name in the result dictionary
            result[variable_name] = os.path.join(root, pattern)
    
    return result

def simulate_pipeline(
                     df_paths,
                     output_dir,
                     simulation_time=250000,
                     gpu_per_node=4,
                     use_mpi=False,
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
        df_paths_rank = df_paths
        device = "0"

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
    output_dir = 'trials/T16_ant3HFM_body4NCO_simulations'
    try:
        os.mkdir(output_dir)
    except:
        pass
    
    # Usage
    output_files = find_fixed_pdb_and_generate_variable('fixed.pdb',
                                             'T16',
                                             directory='trials/T16_ant3HFM_body4NCO')
    
    df = pd.DataFrame(list(output_files.items()),
                       columns = ['VariableName', 'Location'])
    output_files_rf_out =find_fixed_pdb_and_generate_variable('fixed_0.pdb',
                                             'T16',
                                             directory='trials/T16_ant3HFM_body4NCO')
    df_rfout = pd.DataFrame(list(output_files_rf_out.items()),
                       columns = ['VariableName', 'Location'])
    
    
    #output_files_T9_rfout = find_fixed_pdb_and_generate_variable(
    #                            'fixed_0.pdb',
    #                            'T9',
    #                            directory='trials/T9_ant3HFM_body4NCO')
    #
    #df_t9_rfout = pd.DataFrame(list(output_files_T9_rfout.items()),
    #                   columns = ['VariableName', 'Location'])
    #
    #output_files_T10_rfout = find_fixed_pdb_and_generate_variable(
    #                            'fixed_0.pdb',
    #                            'T10',
    #                            directory='trials/T10_ant3HFM_body4NCO')
    #
    #df_t10_rfout = pd.DataFrame(list(output_files_T10_rfout.items()),
    #                   columns = ['VariableName', 'Location'])
    
    df = pd.concat([df, df_rfout])
    df.to_csv('all_pdbs_save_loc.t16.csv', index=False)

if False:
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()
    device = rank % gpu_per_node
    chunks = np.array_split(df, 3)

#cdrlist_rank = np.array_split(cdrlist, int(size))[int(rank)]



if False:
    input_pdb = 'trials/T14_ant3HFM_body4NCO/0/0/1/chaiout/fixed.pdb'
    simulation_time = 250000
    output_dcd = 'output_test_sim_imp/t14_0_0_1_start.dcd'
    output_pdb = 'output_test_sim_imp/t14_0_0_1_start.pdb'
    output_log = 'output_test_sim_imp/t14_0_0_1_start.log'
    d_ind = "0"
    potential_en = simulate_struct(
                        input_pdb,
                        simulation_time,
                        output_dcd,
                        output_pdb,
                        output_log,
                        d_ind=d_ind,
                    )

    print(potential_en)
