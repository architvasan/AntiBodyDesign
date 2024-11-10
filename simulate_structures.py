import os
import pandas as pd
from mpi4py import MPI
import json
import numpy as np
from tqdm import tqdm
import src.simulate.simulation_funcs as sfuncs
import os

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

output_dir = 'trials/T14_ant3HFM_body4NCO_simulations'
try:
    os.mkdir(output_dir)
except:
    pass

# Usage
output_files = find_fixed_pdb_and_generate_variable()
print(output_files)
for variable_name, pdb_path in output_files.items():
    print(f"{variable_name}: {pdb_path}")

df = pd.DataFrame(list(output_files.items()),
                   columns = ['VariableName', 'Location'])

df.to_csv('t14_pdbs_save_loc.csv', index=False)




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