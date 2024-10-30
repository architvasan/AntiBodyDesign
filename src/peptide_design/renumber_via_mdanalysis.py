import MDAnalysis as mda
from MDAnalysis.lib.mdamath import make_whole
from MDAnalysis.lib.util import unique_int_1d

# Load the PDB file
u = mda.Universe('6w4h_cleaned.pdb')

#bad_resnames = ["BGLN", "BGLU", "BLYS", "BASN", "BTHR", "BARG", "BLEU", "BASP", "BVAL", "BSER", "BILE", "BTYR"]
#changefrom_resnames = {"AGLN": "GLN", "AGLU":"GLU", "ALYS":"LYS", "AASN":"ASN", "ATHR":"THR", "AARG":"ARG", "ALEU":"LEU", "AASP":"ASP", "AVAL":"VAL", "ASER":"SER", "AILE":"ILE", "ATYR":"TYR"}
#
#print(u.atoms.select_atoms('resname AGLN'))
#for br in bad_resnames:
#    # Select residues with resname 'BGLN' and remove them
#    print(u.atoms.select_atoms(f'resname {br}'))
#    u.atoms.select_atoms(f'resname {br}').residues.remove()
#
## Rename 'AGLN' to 'GLN'
#for res in u.residues:
#    rnm = res.resname
#    if rnm in changefrom_resnames:#'AGLN':
#        res.resname = changefrom_resnames[rnm]
# Select the residues in chain B

chain_A = u.select_atoms('chainID A and protein and not resname SAM').residues
print(chain_A)
chain_B = u.select_atoms('chainID B and protein and not resname SAM').residues
print(chain_B)

# Renumber residues starting from 1
for i, res in enumerate(chain_A):
    print(res.resid)
    res.resid = i + 1

for i, res in enumerate(chain_B):
    res.resid = i + 1

# Save the modified structure
u.atoms.write('6w4h_renumbered.pdb')

# Write chain A to a separate PDB file
chain_A_sel = u.select_atoms('chainID A and protein and not resname SAM')
chain_A_sel.write('6w4h_chainA.pdb')

# Write chain B to a separate PDB file
chain_B_sel = u.select_atoms('chainID B and protein and not resname SAM')
chain_B_sel.write('6w4h_chainB.pdb')

## Renumber residues starting from 1
#for i, res in enumerate(u.residues):
#    res.resid = i + 1
#
## Save the modified structure
#u.atoms.write('output.pdb')

