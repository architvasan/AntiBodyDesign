import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import json

import ast

if False:
    # Load the file with curly brackets into a dictionary
    with open("dict_t3_r0.dat", "r") as file:
        content = file.read()
        data = ast.literal_eval(content)
    
    #print(data)  # data is now a dictionary
    df = pd.DataFrame(data)
    #print(df)
    
    for r in range(1, 4):
        # Load the file with curly brackets into a dictionary
        with open(f"dict_t3_r{r}.dat", "r") as file:
            content = file.read()
            data = ast.literal_eval(content)
        
        #print(data)  # data is now a dictionary
        df = pd.concat([df, pd.DataFrame(data)])
        
    
    
    for r in range(0,4):
        df = pd.concat([df, pd.read_csv(f'trials/T4/{r}/logout/0/seq_0_heavy_cdr3.csv')])
    df.to_csv('rf_diff_energies.csv', columns=['cdrseq', 'coulen', 'ljen', 'mden'] ,index=False)

df = pd.read_csv('rf_diff_energies.csv')
df['index'] = [it for it in range(len(df))]

df.plot.scatter(x = 'index', y = 'coulen', c='blue')
plt.axhline(y=df['coulen'][0], linestyle='--', label='wt coulombic energy')
plt.savefig('coulomb_energies.png', bbox_inches='tight', dpi=300)
plt.close()

df.plot.scatter(x = 'index', y = 'ljen', c='green')
plt.axhline(y=df['ljen'][0], linestyle='--', label='wt lennard jones energy')
plt.ylim(-100, 20000)
plt.savefig('lenjones_energies.png', bbox_inches='tight', dpi=300)
plt.close()

df.plot.scatter(x = 'index', y = 'mden', c='black')
plt.axhline(y=df['mden'][0], linestyle='--', label='wt md energy')
plt.savefig('md_energies.png', bbox_inches='tight', dpi=300)
plt.close()


# Load the file with curly brackets into a dictionary
#with open("dict_t3_r0.dat", "r") as file:
#    data = json.load(file)
#
#print(data)  # data is now a dictionary

#with open("dict_t3_r0.dat", "r") as fil:
#    dict_0_dat = fil.read()
#    print(dict(dict_0_dat[0]))
#    dict_df = pd.DataFrame(dict_0_dat) 
#
#print(dict_df)
