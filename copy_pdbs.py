import os
import shutil

# Define the input file and output directory
#input_file = 'all_pdbs_save_loc.t16.csv'
#output_dir = 'sim_pdbs'
input_file = 'all_pdbs_save_loc_t2_afterdpo.csv'
output_dir = 'sim_pdbs'

# Ensure the output directory exists
os.makedirs(output_dir, exist_ok=True)

# Open the input file and process each line
with open(input_file, 'r') as f:
    for line in f:
        # Split the line into two parts: name and path
        name, path = line.strip().split(',')
        
        # Define the new file name and path in the output directory
        output_path = os.path.join(output_dir, f"{name}.fixed.pdb")
        
        # Copy the file to the output directory with the new name
        try:
            shutil.copy(path, output_path)
            print(f"Copied {path} to {output_path}")
        except FileNotFoundError:
            print(f"File {path} not found. Skipping.")
        except Exception as e:
            print(f"Error copying {path} to {output_path}: {e}")
