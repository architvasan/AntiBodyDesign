import pandas as pd
import argparse

def compare_column(file1, file2, column_to_compare):
    """
    Compare a specific column between two CSV files and find differing rows.

    Args:
        file1 (str): Path to the first CSV file.
        file2 (str): Path to the second CSV file.
        column_to_compare (str): The column to compare between the two files.

    Returns:
        pd.DataFrame: DataFrame containing rows from file1 where the column differs.
    """
    # Load the two DataFrames
    df1 = pd.read_csv(file1)
    df2 = pd.read_csv(file2)
    df2 = df2.rename(columns={'cdrseq': 'generated_seqs'} )
    df2 = df2.reindex(df1.index)
    # Ensure the column exists in both DataFrames
    if column_to_compare not in df1.columns or column_to_compare not in df2.columns:
        raise ValueError(f"Column '{column_to_compare}' not found in one or both files.")

    # Compare the specified column
    differences = df1[df1[column_to_compare] != df2[column_to_compare]]

    return differences

for i in range(8):
    orig_loc = 'iedb_tables_pos_and_neg'
    trial_loc = 'trials/T1_ant3HFM_body4NCO_afterdpo_norf'
    diff_df = compare_column(f'{orig_loc}/generated_seqs_dpo_{i}.csv', f'{trial_loc}/{i}/0/diff_cdr_results.csv', 'generated_seqs')
    print(diff_df)
    diff_df.to_csv(f'{orig_loc}/generated_seqs_dpo_{i}_rst0.csv', columns = ['generated_seqs', 'probs_ratio'], index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compare a column between two CSV files and find differences.")
    parser.add_argument("--file1", type=str, help="Path to the first CSV file.")
    parser.add_argument("--file2", type=str, help="Path to the second CSV file.")
    parser.add_argument("--column", type=str, help="Name of the column to compare.")

    args = parser.parse_args()

    try:
        # Run the comparison
        differing_rows = compare_column(args.file1, args.file2, args.column)
        
        # Output the result
        if differing_rows.empty:
            print("No differences found in the specified column.")
        else:
            print("Differences found in the following rows:")
            print(differing_rows)
    except Exception as e:
        print(f"Error: {e}")

