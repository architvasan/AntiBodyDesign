import pandas as pd
import numpy as np

df = pd.read_csv('trials/T14_ant3HFM_body4NCO/0/df_cdrs_energies.csv')
for i in range(1,40):
    df = pd.concat([df, pd.read_csv(f'trials/T14_ant3HFM_body4NCO/{i}/df_cdrs_energies.csv')])

df.to_csv(f't14_all_cdrs_ant3hfm_body4nco.csv', index=False)

df = pd.read_csv('trials/T10_ant3HFM_body4NCO/0/0/diff_cdr_results.csv')
for i in range(1,8):
    df = pd.concat([df, pd.read_csv(f'trials/T14_ant3HFM_body4NCO/{i}/0/diff_cdr_results.csv')])

df.to_csv(f't10_all_cdrs_ant3hfm_body4nco.csv', index=False)

df = pd.read_csv('trials/T9_ant3HFM_body4NCO/0/0/diff_cdr_results.csv')
for i in range(1,8):
    df = pd.concat([df, pd.read_csv(f'trials/T9_ant3HFM_body4NCO/{i}/0/diff_cdr_results.csv')])

df.to_csv(f't9_all_cdrs_ant3hfm_body4nco.csv', index=False)