import argparse
from Bio.PDB import PDBList

def download_pdb(pdb_id, save_dir="."):
    # Create a PDBList object
    pdb = PDBList()
    # Download the PDB file
    pdb_file = pdb.retrieve_pdb_file(pdb_id, file_format="pdb", pdir=save_dir)
    print(f"PDB file downloaded to: {pdb_file}")

if __name__ == "__main__":
    # Set up argument parsing
    parser = argparse.ArgumentParser(description="Download a PDB file using its PDB ID.")
    parser.add_argument("pdb_id", type=str, help="The 4-character PDB ID of the structure to download.")
    parser.add_argument(
        "--save_dir", type=str, default=".", help="Directory to save the PDB file (default: current directory)."
    )

    # Parse the arguments
    args = parser.parse_args()

    # Call the download function with the provided arguments
    download_pdb(args.pdb_id, args.save_dir)

