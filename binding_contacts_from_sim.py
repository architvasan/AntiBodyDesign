import src.utils.seq_frompdb as seq_frompdb
import pandas as pd
import numpy as np
import glob
import src.analyze.contacts_per_frame as contacts_analysis
from tqdm import tqdm
import MDAnalysis as mda
from tempfile import NamedTemporaryFile
import src.analyze.lie_2parts as lie
def find_substrings(sequence, substrings):
    """Find which substrings are present in the given sequence."""
    return [s for s in substrings if s in sequence]

def _obtain_cdr_pdb(pdb, inp_data_file):
    df = pd.read_csv(inp_data_file)
    df_cdrs = list(df['cdrseq'])
    print(df_cdrs)
    seq_pdb = seq_frompdb.get_seq_from_pdb(pdb)
    print(seq_pdb)
    cdrs = find_substrings(seq_pdb[0], df_cdrs)
    return cdrs

def obtain_cdr_pdb(pdb, cdrlist):
    seq_pdb = seq_frompdb.get_seq_from_pdb(pdb)
    cdrs = find_substrings(seq_pdb[0], cdrlist)
    return cdrs, seq_pdb

def get_dcd_files(directory):
    return glob.glob(f"{directory}/**/*.dcd", recursive=True)

def parse_dir(dcd_directory, pdb_directory):
    all_dcds = get_dcd_files(dcd_directory)
    all_pdbs = [f'{pdb_directory}/{s.rsplit(".dcd")[0].rsplit("/")[-1]}.fixed.pdb' for s in all_dcds]
    all_logs = [f'{s.rsplit(".dcd")[0]}.log' for s in all_dcds]
    #print(all_dcds)
    #print(all_pdbs)
    return all_pdbs, all_dcds, all_logs

def lie_on_all_frames(
                      pdb_file,
                      dcd_file):

    # Load the trajectory in MDAnalysis
    universe = mda.Universe(pdb_file, dcd_file)
    coul_en_list = []
    lj_en_list = []
    for frame in (universe.trajectory[-1:]):
        # Create a temporary PDB file to hold the current frame coordinates
        with NamedTemporaryFile(suffix=".pdb") as temp_pdb:
            temp_pdb_filename = temp_pdb.name
            universe.atoms.write(temp_pdb_filename)
            #print(f"Temporary PDB saved for frame {frame.frame}: {temp_pdb_filename}")

            coul_en, lj_en = lie.lie(temp_pdb_filename, 'A', 1, 500, 'C')
            coul_en_list.append(coul_en)
            lj_en_list.append(lj_en)
    del(universe)
    return coul_en_list, lj_en_list
def parse_cdrs_dir(dcd_directory, pdb_directory, inp_data_file):
    all_pdbs, all_dcds, all_logs = parse_dir(dcd_directory, pdb_directory)
    df = pd.read_csv(inp_data_file)
    cdrlist = list(df['cdrseq'])
    kcal_kj_conv = 0.239005736
    dict_cdr = {}
    dict_cdr['cdrseq'] = []
    dict_cdr['mden'] = []
    dict_cdr['contacts'] = []
    dict_cdr['coul_int_en'] = []
    dict_cdr['lj_int_en'] = []
    for pdbit, dcdit, logit in tqdm(zip(all_pdbs, all_dcds, all_logs)):
        try:
            seq_cdr, seqpdb = obtain_cdr_pdb(pdbit, cdrlist) 
            rid_init = "".join(seqpdb).find(seq_cdr[0])
            rid_fin = rid_init + len(seq_cdr[0])
            log_df = pd.read_csv(logit)
            av_en = kcal_kj_conv * np.mean(list(log_df['#"Potential Energy (kJ/mole)"']))
            conts_per_frame = contacts_analysis.main(pdbit, 
                                dcdit,
                                f'segid A and name CA',
                                f'segid C and name CA',
                                8,
                                'out_contacts_per_frame',
                                )
            coul_perf, lj_perf = lie_on_all_frames(pdbit, dcdit)
            print(coul_perf)
            print(lj_perf)
            # and resid {rid_init}-{rid_fin} 
            dict_cdr['cdrseq'].append(seq_cdr[0])
            dict_cdr['mden'].append(av_en)
            dict_cdr['contacts'].append(np.mean(conts_per_frame))
            dict_cdr['coul_int_en'].append(np.mean(coul_perf))
            dict_cdr['lj_int_en'].append(np.mean(lj_perf))

            #print(seq_cdr)
            #print(av_en) 

        except Exception as e:
            print(e)
            continue
    return dict_cdr


c_en_init, l_en_init = lie_on_all_frames(pdb_file='3hfm_outs/fixed_0.pdb', dcd_file = '3hfm_outs/3hfm.dcd')
print(c_en_init)
print(l_en_init)
# Example usage
dcd_directory = "trials/T16_ant3HFM_body4NCO_simulations/"
pdb_directory = 'sim_pdbs'

dict_cdr_t16 = parse_cdrs_dir(dcd_directory, 
               pdb_directory, 
               't16_all_cdrs_ant3hfm_body4nco.csv'
               )

df_cdr_t16 = pd.DataFrame(dict_cdr_t16)
df_cdr_t16.to_csv('t16_gen_cdrs_mden.csv', index=False)

dcd_directory = "trials/T9_10_14_antdy4NCO_simulations"
pdb_directory = 'sim_pdbs'

dict_cdr_t9 = parse_cdrs_dir(dcd_directory, 
               pdb_directory, 
               't9_all_cdrs_ant3hfm_body4nco.csv'
               )
df_cdr_t9 = pd.DataFrame(dict_cdr_t9)
df_cdr_t9.to_csv('t9_gen_cdrs_mden.csv', index=False)

dict_cdr_t10 = parse_cdrs_dir(dcd_directory, 
               pdb_directory, 
               't10_all_cdrs_ant3hfm_body4nco.csv'
               )
df_cdr_t10 = pd.DataFrame(dict_cdr_t10)
df_cdr_t10.to_csv('t10_gen_cdrs_mden.csv', index=False)

dict_cdr_t14 = parse_cdrs_dir(dcd_directory, 
               pdb_directory, 
               't14_all_cdrs_ant3hfm_body4nco.csv'
               )
df_cdr_t14 = pd.DataFrame(dict_cdr_t14)
df_cdr_t14.to_csv('t14_gen_cdrs_mden.csv', index=False)

df_cdr_all = pd.concat([df_cdr_t9, df_cdr_t10, df_cdr_t14, df_cdr_t16])
df_cdr_all.to_csv('all_gen_cdrs_mden.csv', index=False)

#dcd_files = get_dcd_files(directory)
#print("DCD files:", dcd_files)




#print(obtain_cdr_pdb('sim_pdbs/T16_1_0_4_rfout_chai_struct_17_output.fixed.pdb', 
                     #'t16_all_cdrs_ant3hfm_body4nco.csv'))
