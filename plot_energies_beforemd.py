import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import json
import ast
import math
import seaborn as sns
def boltzmann_term(
                    en,
                    en_mean,
                    en_std,
                    kbT=0.592,
                    ):
    #print(en)
    #print(en_mean)
    #print(en_std)
    #print(-((en-en_mean))/(en_std))
    try:
        #print(-en)
        #print(math.exp(1))
        #print(math.exp(-en/kbT))
        #return math.exp(-(en)/kbT)
        return math.exp(-((en-en_mean))/(en_std))
    except:
        return 1

### Before dpo
df = pd.read_csv('trials/T16_ant3HFM_body4NCO/0/0/diff_cdr_results.csv')
for i in range(1, 8):
    df = pd.concat([df, pd.read_csv(f'trials/T16_ant3HFM_body4NCO/{i}/0/diff_cdr_results.csv')])

df = df[df['coulen']!=0]
df = df[df['ljen']<1000]
### After dpo
df_after = pd.read_csv('trials/T2_ant3HFM_body4NCO_afterdpo_norf/0/0/diff_cdr_results.csv')
for i in range(1, 8):
    df_after = pd.concat([df_after, pd.read_csv(f'trials/T2_ant3HFM_body4NCO_afterdpo_norf/{i}/0/diff_cdr_results.csv')])

print(df_after)
print(df)
df = df.reset_index()
df_after = df_after.reset_index()
df_after = df_after[df_after['coulen']!=0]
df_after = df_after[df_after['ljen']<1000]

#ax1 = df.plot.scatter(x = 'index', y = 'coulen', c='blue')
#df_after.plot.scatter(x = 'index', y = 'coulen', c='red', ax = ax1)
sns.kdeplot(df['coulen'], color='b', shade=True)
sns.kdeplot(df_after['coulen'], color='r', shade=True)
plt.legend()
plt.savefig('coulen_beforemd_dpocomp.png', bbox_inches='tight', dpi=300)
plt.close()

# ax2 = df.plot.scatter(x = 'index', y = 'ljen', c='blue')
# df_after.plot.scatter(x = 'index', y = 'ljen', c='red', ax = ax2)
sns.kdeplot(df['ljen'], color='b', shade=True)
sns.kdeplot(df_after['ljen'], color='r', shade=True)
plt.legend()
plt.savefig('ljen_beforemd_dpocomp.png', bbox_inches='tight', dpi=300)
plt.close()

# ax3 = df.plot.scatter(x = 'index', y = 'mden', c='blue')
# df_after.plot.scatter(x = 'index', y = 'mden', c='red', ax=ax3)

sns.kdeplot(df['mden'], color='b', shade=True)
sns.kdeplot(df_after['mden'], color='r', shade=True)

plt.legend()
plt.savefig('mden_beforemd_dpocomp.png', bbox_inches='tight', dpi=300)
plt.close()
