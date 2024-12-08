import numpy as np
import json
import pandas as pd

df = pd.read_csv('./all_pdbs_save_loc_t2_afterdpo.csv')
df_split = np.array_split(df, 8)
for i in range(0,len(df_split)):
    df_split[i].to_csv(f'all_pdbs_save_loc_t2_afterdpo_{i}.csv', index=False)

if False:
    with open('./all_pdbs_save_loc.t16.csv','r') as file:
        cdrlist = json.load(file)
    
    cdrlist_set = set(cdrlist['all_train_sequences'])
    cdrlist_uniq = list(cdrlist_set)
    with open(f'iedb_tables_pos_and_neg/heavy-cdr3_all.json', 'w') as file:
        json.dump({'all_train_sequences': list(cdrlist_uniq)}, file)
    
    cdrlist_split = np.array_split(cdrlist_uniq, 8)
    for i in range(8):
        with open(f'iedb_tables_pos_and_neg/heavy-cdr3_rank{i}.json', 'w') as file:
            json.dump({'all_train_sequences': list(cdrlist_split[i])}, file)
