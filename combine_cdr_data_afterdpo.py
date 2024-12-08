import pandas as pd
import numpy as np

df = pd.read_csv('trials/T2_ant3HFM_body4NCO_afterdpo_norf/0/0/diff_cdr_results.csv')
for i in range(1,8):
    df = pd.concat([df, pd.read_csv(f'trials/T2_ant3HFM_body4NCO_afterdpo_norf/{i}/0/diff_cdr_results.csv')])

df.to_csv(f't2_all_cdrs_ant3hfm_body4nco_afterdpo.csv', index=False)
