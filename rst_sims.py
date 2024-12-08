import os
import pandas as pd

def rst_csv(csv_file_path,
            search_column,
            directory_path,
            output_csv_path):

    # Parameters
    """
    csv_file_path = "your_data.csv"  # Path to your CSV file
    search_column = "pattern_column"  # Column in the DataFrame to match against filenames
    directory_path = "your_directory"  # Directory to search for matching patterns
    output_csv_path = "filtered_data.csv"  # Path to save the filtered DataFrame
    """

    # Step 1: Load the CSV file into a Pandas DataFrame
    df = pd.read_csv(csv_file_path)

    # Step 2: Get the list of base filenames (without `.pdb`) in the specified directory
    filenames = set(
                os.path.splitext(filename)[0] for filename in os.listdir(directory_path) if filename.endswith(".pdb")
                )

    # Step 3: Remove rows from the DataFrame where the column matches any filename
    df_filtered = df[~df[search_column].isin(filenames)]

    # Step 4: Save the filtered DataFrame to a new CSV file
    df_filtered.to_csv(output_csv_path, index=False)

    print(f"Filtered DataFrame saved to {output_csv_path}")
    return

for i in range(0,8):
    csv_file_path = f'all_pdbs_save_loc_t2_afterdpo_{i}.csv'
    search_column = 'VariableName'
    directory_path = f'trials/T2_ant3HFM_body4NCO_simulations_afterdpo/{i}'
    output_csv_path = f'all_pdbs_save_loc_t2_afterdpo_{i}.rst0.csv'
    rst_csv(csv_file_path,
            search_column,
            directory_path,
            output_csv_path)