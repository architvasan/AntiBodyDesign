import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import json
import ast
import math

def boltzmann_term(
                    en,
                    en_mean,
                    en_std,
                    kbT=0.592,
                    ):

    try:
        return math.exp(-((en-en_mean))/(en_std))
    except:
        return 1

df = pd.read_csv('all_gen_cdrs_mden_aftermd.csv')
df = df.drop_duplicates(subset='cdrseq')
df['cen_ljen'] = df['coul_int_en'] + df['lj_int_en']

df_dpo = pd.read_csv('t2_afterdpo_gen_cdrs_mden.csv')
df_dpo = df_dpo[df_dpo['coul_int_en']!=0]
#df_dpo = df_dpo.drop_duplicates(subset='cdrseq')
df_dpo['cen_ljen'] = df_dpo['coul_int_en'] + df_dpo['lj_int_en']

df_3hfm = pd.read_csv('initial_3hfm_mden_lie_contacts.dat')
inten_3hfm = df_3hfm['ljen'][0] + df_3hfm['coulen'][0]
df_3hfm['cen_ljen'] = inten_3hfm

df_en_tot = list(df['cen_ljen'])
print(len(df_en_tot))
df_en_tot.extend(list(df_dpo['cen_ljen']))
print(len(df_en_tot))
df_en_tot.append(inten_3hfm)
print(len(df_en_tot))
int_en_mean = np.mean(df_en_tot)
int_en_std = np.std(df_en_tot)

boltz_terms = [boltzmann_term(list(df['cen_ljen'])[it],int_en_mean, int_en_std) for it in range(len(df))]
boltz_term_3hfm = boltzmann_term(inten_3hfm, int_en_mean, int_en_std)
boltz_terms_dpo = [boltzmann_term(list(df_dpo['cen_ljen'])[it],int_en_mean, int_en_std) for it in range(len(df_dpo))]

boltz_terms_main = boltz_terms.copy()
boltz_terms_main.extend(boltz_terms_dpo)
boltz_terms_main.append(boltz_term_3hfm)
part = np.sum(boltz_terms_main)

df['enweight'] = boltz_terms/part
df_dpo['enweight'] = boltz_terms_dpo/part
enweight_3hfm = boltz_term_3hfm/part
#df_3hfm.to_csv('3hfm_cdr_energy.csv', index=False)

contacts_3hfm = list(df_3hfm['contacts'])[0]
en_3hfm = list(df_3hfm['mden'])[0]
coulen_3hfm = list(df_3hfm['coulen'])[0]
ljen_3hfm = list(df_3hfm['ljen'])[0]

df['contacts_weighted'] = part * df['contacts'] * df['enweight']
print(min(df_dpo['enweight']))
df_dpo['contacts_weighted'] = part * df_dpo['contacts'] * df_dpo['enweight']
df_dpo = df_dpo[df_dpo['contacts_weighted']>0]
df_3hfm_contacts_weighted = part * df_3hfm['contacts'][0]*enweight_3hfm
df_3hfm_conts = df_3hfm['contacts']

df_dpo.to_csv('all_cdr_afterdpo_energies.csv', index=False)

print(df_dpo[df_dpo['contacts_weighted']>df_3hfm_contacts_weighted])
df_dpo_pos = df_dpo[df_dpo['contacts_weighted']>df_3hfm_contacts_weighted]
print(f"Percentage positive: {len(df_dpo_pos)/len(df_dpo)}")
import seaborn as sns
sns.kdeplot(df['coul_int_en'], color='b', shade=True, alpha=0.6)
sns.kdeplot(df_dpo['coul_int_en'], color='r', shade=True, alpha=0.6)
plt.axvline(x=coulen_3hfm, linestyle='--', color = 'black', label='wt/3hfm coulombic energy')
plt.legend()
plt.savefig('coulen_aftermd_dpocomp.png', bbox_inches='tight', dpi=300)
plt.close()

sns.kdeplot(df['lj_int_en'], color='b', shade=True, alpha=0.6)
sns.kdeplot(df_dpo['lj_int_en'], color='r', shade=True, alpha=0.6)
plt.axvline(x=ljen_3hfm, linestyle='--', color = 'black', label='wt/3hfm lennard jones energy')
plt.legend()
plt.savefig('ljen_aftermd_dpocomp.png', bbox_inches='tight', dpi=300)
plt.close()

sns.histplot(df['contacts_weighted'], color='b', bins=50, stat='probability')#, shade=True)
sns.histplot(df_dpo['contacts_weighted'], color='r', bins=1000, stat='probability')#, shade=True)
plt.axvline(x=df_3hfm_contacts_weighted, linestyle='--', color = 'black', label='wt/3hfm weighted contacts')
plt.legend()
plt.xlim(0,1000)
plt.savefig('contacts_weighted_aftermd_dpocomp.png', bbox_inches='tight', dpi=300)
plt.close()

sns.histplot(df['contacts_weighted'], color='b', bins=5, stat='probability', alpha=0.8)#, shade=True)
sns.histplot(df_dpo['contacts_weighted'], bins=5, color='r', stat='probability', alpha=0.8)#, shade=True)
plt.axvline(x=df_3hfm_contacts_weighted, linestyle='--', color = 'black', label='wt/3hfm weighted contacts')
plt.legend()
plt.ylim(0,0.05)
plt.xlim(df_3hfm_contacts_weighted,1000)
plt.savefig('contacts_weighted_aftermd_dpocomp_zoomin.png', bbox_inches='tight', dpi=300)
plt.close()

sns.histplot(df['contacts'], color='b', stat='probability')#, shade=True)
sns.histplot(df_dpo['contacts'], color='r', stat='probability')#, shade=True)
#plt.axvline(x=df_3hfm_conts[0], linestyle='--', color = 'black', label='wt/3hfm weighted contacts')
plt.legend()
plt.xlim(0,100)
plt.savefig('contacts_aftermd_dpocomp.png', bbox_inches='tight', dpi=300)
plt.close()

if False:
    # df['index'] = [it for it in range(len(df))]
    # print(df[df['contacts_weighted']>contacts_weighted_3hfm])
    # df_pos = df[df['contacts_weighted']>contacts_weighted_3hfm]
    # df_neg = df[df['contacts_weighted']<=contacts_weighted_3hfm]
    
    # df_pos.to_csv('cdr_gen_ens_pos.csv', index=False)
    # df_neg.to_csv('cdr_gen_ens_neg.csv', index=False)
    
    ######## After DPO ####################
    
    
    df_3hfm = pd.read_csv('initial_3hfm_mden_lie_contacts.dat')
    inten_3hfm = df_3hfm['ljen'][0] + df_3hfm['coulen'][0]
    df_3hfm['cen_ljen'] = inten_3hfm
    df_en_tot_dpo = list(df_dpo['cen_ljen']).append(inten_3hfm)
    int_en_mean_dpo = np.mean([np.mean(df_dpo['cen_ljen']), np.mean(df['cen_ljen'])])
    int_en_std_dpo = np.std(np.concatenate((df_dpo['cen_ljen'], df['cen_ljen'])))
    df_dpo = df_dpo.drop_duplicates(subset='cdrseq')
    #df.to_csv('all_cdr_energies.csv', index=False)
    
    boltz_terms_dpo = [boltzmann_term(list(df_dpo['cen_ljen'])[it],int_en_mean_dpo, int_en_std_dpo) for it in range(len(df_dpo))]
    boltz_term_3hfm = boltzmann_term(inten_3hfm, int_en_mean_dpo, int_en_std_dpo)
    boltz_terms_main_dpo = boltz_terms_dpo.copy()
    boltz_terms_main_dpo.append(boltz_term_3hfm)
    part = np.sum(boltz_terms_main) + np.sum(boltz_terms_main_dpo)
    
    df['enweight'] = boltz_terms/part
    enweight_3hfm = boltz_term_3hfm/part
    
    #print("enweights + contacts are:")
    #print(df['enweight'])
    #print(df['contacts'])
    #print(df['contacts']*df['enweight'])
    df['contacts_weighted'] = part * df['contacts'] * df['enweight'] #[list(df['enweight'])[it] * list(df['contacts'])[it] for it in range(len(df))]
    
    contacts_weighted_3hfm = part * df_3hfm['contacts'][0]*enweight_3hfm
    df_3hfm['contacts_weighted'] = contacts_weighted_3hfm
    print(f"weighted conts for 3hfm: {contacts_weighted_3hfm}")
    #df_3hfm['contacts_weighted'] = df_3hfm['contacts']#boltz_terms/part
    
    
    
    df_dpo['enweight'] = boltz_terms_dpo/part_dpo
    enweight_3hfm = boltz_term_3hfm/part_dpo
    
    #print("enweights + contacts are:")
    #print(df['enweight'])
    #print(df['contacts'])
    #print(df['contacts']*df['enweight'])
    df_dpo['contacts_weighted'] = part_dpo * df_dpo['contacts'] * df_dpo['enweight'] #[list(df['enweight'])[it] * list(df['contacts'])[it] for it in range(len(df))]
    
    contacts_weighted_3hfm = part_dpo * df_3hfm['contacts'][0]*enweight_3hfm
    df_3hfm['contacts_weighted'] = contacts_weighted_3hfm
    print(f"weighted conts for 3hfm: {contacts_weighted_3hfm}")
    #df_3hfm['contacts_weighted'] = df_3hfm['contacts']#boltz_terms/part
    
    df_dpo = df_dpo.drop_duplicates(subset='cdrseq')
    df_dpo.to_csv('all_cdr_afterdpo_energies.csv', index=False)
    
    df_3hfm.to_csv('3hfm_cdr_energy.csv', index=False)
    contacts_3hfm = list(df_3hfm['contacts'])[0]
    en_3hfm = list(df_3hfm['mden'])[0]
    coulen_3hfm = list(df_3hfm['coulen'])[0]
    ljen_3hfm = list(df_3hfm['ljen'])[0]
    
    df['index'] = [it for it in range(len(df))]
    print(df[df['contacts_weighted']>contacts_weighted_3hfm])
    df_pos = df[df['contacts_weighted']>contacts_weighted_3hfm]
    df_neg = df[df['contacts_weighted']<=contacts_weighted_3hfm]
    
    df_pos.to_csv('cdr_gen_ens_pos.csv', index=False)
    df_neg.to_csv('cdr_gen_ens_neg.csv', index=False)
    
    
    
    
    df.plot.scatter(x = 'index', y = 'contacts_weighted', c='orange')
    plt.axhline(y=contacts_weighted_3hfm, linestyle='--', color='red', label='wt/3hfm weighted contacts')
    #plt.ylim(0,100)
    plt.legend()
    plt.savefig('contact_enweighted.lysozyme.png', bbox_inches='tight', dpi=300)
    plt.close()
    
    df.plot.scatter(x = 'index', y = 'coul_int_en', c='blue')
    plt.axhline(y=coulen_3hfm, linestyle='--', color='red', label='wt/3hfm coulombic energy')
    plt.legend()
    plt.savefig('coulomb_energies.lysozyme.png', bbox_inches='tight', dpi=300)
    plt.close()
    
    df.plot.scatter(x = 'index', y = 'lj_int_en', c='green')
    plt.axhline(y=ljen_3hfm, linestyle='--', color = 'red', label='wt/3hfm lennard jones energy')
    #plt.ylim(-100, 200)
    plt.legend()
    plt.savefig('lenjones_energies.lysozyme.png', bbox_inches='tight', dpi=300)
    plt.close()
    
    df.plot.scatter(x = 'index', y = 'mden', c='black')
    plt.axhline(y=en_3hfm, linestyle='--', color='red', label='wt/3hfm md energy')
    plt.legend()
    plt.savefig('md_energies.lysozyme.png', bbox_inches='tight', dpi=300)
    plt.close()
    
    df.plot.scatter(x = 'index', y = 'contacts', c='purple')
    plt.axhline(y=contacts_3hfm, linestyle='--', color = 'red', label='wt/3hfm contacts')
    plt.legend()
    plt.savefig('contacts.lysozyme.png', bbox_inches='tight', dpi=300)
    plt.close()