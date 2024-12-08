import os
import pandas as pd
import json
import numpy as np
from tqdm import tqdm
import os

def find_fixed_pdb_and_generate_variable(pattern,
                                         outpatt,
                                         directory):
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

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-pc',
                        '--pattern_chai',
                        type=str,
                        help='filename pattern to search for initial chai structs')

    parser.add_argument('-op',
                        '--outpatt',
                        type=str,
                        help='filename pattern to save in df')
    
    parser.add_argument('-D',
                        '--directory',
                        type=str,
                        help='directory to search in')

    parser.add_argument('-of',
                        '--outfile',
                        type=str,
                        help='csvfile to save pdb locations')
    
    args = parser.parse_args()

    output_files = find_fixed_pdb_and_generate_variable(
                                args.pattern_chai,
                                args.outpatt,
                                directory=args.directory)

    df = pd.DataFrame(list(output_files.items()),
                   columns = ['VariableName', 'Location'])
    
    df.to_csv(args.outfile, index=False)



    
    # output_files = find_fixed_pdb_and_generate_variable('fixed.pdb',
    #                                          'T16',
    #                                          directory='trials/T16_ant3HFM_body4NCO')
    
    # df = pd.DataFrame(list(output_files.items()),
    #                    columns = ['VariableName', 'Location'])
    # output_files_rf_out =find_fixed_pdb_and_generate_variable('fixed_0.pdb',
    #                                          'T16',
    #                                          directory='trials/T16_ant3HFM_body4NCO')
    # df_rfout = pd.DataFrame(list(output_files_rf_out.items()),
    #                    columns = ['VariableName', 'Location'])
    


# if True:
#     output_dir = 'trials/T16_ant3HFM_body4NCO_simulations'
#     try:
#         os.mkdir(output_dir)
#     except:
#         pass
    
#     # Usage
#     output_files = find_fixed_pdb_and_generate_variable('fixed.pdb',
#                                              'T16',
#                                              directory='trials/T16_ant3HFM_body4NCO')
    
#     df = pd.DataFrame(list(output_files.items()),
#                        columns = ['VariableName', 'Location'])
#     output_files_rf_out =find_fixed_pdb_and_generate_variable('fixed_0.pdb',
#                                              'T16',
#                                              directory='trials/T16_ant3HFM_body4NCO')
#     df_rfout = pd.DataFrame(list(output_files_rf_out.items()),
#                        columns = ['VariableName', 'Location'])
    
    
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
    
    # df = pd.concat([df, df_rfout])
    # df.to_csv('all_pdbs_save_loc.t16.csv', index=False)

# if False:
#     comm = MPI.COMM_WORLD
#     size = comm.Get_size()
#     rank = comm.Get_rank()
#     device = rank % gpu_per_node
#     chunks = np.array_split(df, 3)

# #cdrlist_rank = np.array_split(cdrlist, int(size))[int(rank)]



# if False:
#     input_pdb = 'trials/T14_ant3HFM_body4NCO/0/0/1/chaiout/fixed.pdb'
#     simulation_time = 250000
#     output_dcd = 'output_test_sim_imp/t14_0_0_1_start.dcd'
#     output_pdb = 'output_test_sim_imp/t14_0_0_1_start.pdb'
#     output_log = 'output_test_sim_imp/t14_0_0_1_start.log'
#     d_ind = "0"
#     potential_en = simulate_struct(
#                         input_pdb,
#                         simulation_time,
#                         output_dcd,
#                         output_pdb,
#                         output_log,
#                         d_ind=d_ind,
#                     )

#     print(potential_en)
