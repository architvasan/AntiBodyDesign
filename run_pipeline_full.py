#import src.utils.id_cdrloop as cdrloop
import src.fold_ai.chai_pred as chai_pred
import pandas as pd
import os
import src.utils.cif2pdb as cif2pdb
import MDAnalysis as mda
import src.utils.fix_rfdiff_out as fix_rfdiff
import src.utils.seq_frompdb as seq_frompdb
import glob
import src.analyze.lie as lie
import src.rfdiffrun.partialdiff_loop as partial_loop
import src.rfdiffrun.fixrfpdb_seq as fixrfpdb_seq
import src.rfdiffrun.silenttools as silenttools
import src.rfdiffrun.dlbinder as dlbinder
import src.analyze.md_energy as md_energy

'''
Steps:
0. Load in pandas dataframe
1. Determine sequencs for light + heavy chains + antigens
2. Determine resid numbers for CDR loop to target
3. Load light + heavy chains into CHAI-1
4. Load CHAI-1 structure into rfdiffusion to generate new CDR loops with partial diff
4.1 Fix residues with known sequences.
4.2 Convert pdbs in directory to .silent file
5 Load rfdiffusion .silent file into dl binder design
6. Load dlbind design sequence into CHAI-1
6.1 Convert .silent files back to .pdb
8. Run a short MD simulation for each structure
9. Determine interaction energy + md energy from each simulation
'''

def run_pipeline_single(rowit,
                        data,
                        cdrloop,
                        chaintarget,
                        chai_dir,
                        rf_dir,
                        rf_script_path,
                        num_designs,
                        se3_env,
                        dlbind_env):
    '''
    step 1
    '''
    data_it = data.loc[rowit]
    data_it_heavy = data_it['heavy_chain']
    data_it_light = data_it['light_chain']
    data_it_ant = data_it['antigen']

    '''
    step 2
    '''
    #data_it_ant = data_it['antigen']
    data_it_target = data_it[chaintarget]
    data_it_patt = data_it[cdrloop]
    rid_init = data_it_target.find(data_it_patt)
    rid_fin = rid_init + len(data_it_patt)

    if False:
        '''
        step 3
        '''
        chai_pred.fold_chai_body_ant(
                            data_it_heavy,
                            data_it_light,
                            data_it_ant,
                            f'{chai_dir}/temp.fasta',
                            chai_dir,
                            device=0
                            )
        '''
        step 3.1
        use cif2pdb util to convert cif to pdb
        '''
        cif2pdb.cif2pdb(f'{chai_dir}/pred.model_idx_0.cif')

    '''
    step 4
    Construct and execute the command
    * Assume for now that we are dealing with heavy chain.
    Will modify for light chain later
    '''

    len_heavy = len(list(data_it_heavy))
    len_light = len(list(data_it_light))
    len_ant = len(list(data_it_ant))
    
    partial_loop.rfdiff_full(
          f'{chai_dir}/pred.model_idx_0.pdb',
          rf_dir,
          rid_init,
          rid_fin,
          len_heavy,
          len_light,
          len_ant,
          rf_script_path,
          se3_env,
          partial_steps=10,
          num_designs=num_designs,
          )

    '''
    step 4.15
    move any residue with > len(heavy_chain) to chain B + > len(heavy_chain+light_chain) to chain C
    '''
    for file in os.listdir(rf_dir):
        print(file)
        if file.endswith(".pdb"):
            print(file)
            pdb_path = os.path.join(rf_dir, file)
            fix_rfdiff.modify_chain_full(pdb_path,
                                         len_heavy,
                                         len_light,
                                        )

    '''
    step 4.1
    fix residues in pdb 
    '''

    fixrfpdb_seq.fixpdb(rf_dir,
                        dlbind_env,
                        '/eagle/projects/datascience/avasan/RFDiffusionProject/dl_binder_design/helper_scripts',
                        )
    
    '''
    step 4.2
    convert pdb to silent
    '''

    silenttools.pdb2silent(dlbind_env,
                            rf_dir,
                            '/eagle/projects/datascience/avasan/RFDiffusionProject/dl_binder_design/include/silent_tools',
                            )


    '''
    step 5
    '''
    dlbinder.protein_mpnn(dlbind_env,
                              '/eagle/projects/datascience/avasan/RFDiffusionProject/dl_binder_design/mpnn_fr',
                              rf_dir)

    '''
    step 6 alternative
    use chai-1 to obtain new structure
    1. convert mpnnout.silent back to pdbs
    2. for each pdb: determine sequence using seq tool
    3. for each sequence: plug sequence into chai-1 to determine final fold
    '''

    '''
    step 6.1
    '''
    silenttools.extractpdb('/eagle/datascience/avasan/Simulations/Antibody_Design',
                           dlbind_env,
                           rf_dir,
                           '/eagle/projects/datascience/avasan/RFDiffusionProject/dl_binder_design/include/silent_tools',
                           )
    
    '''
    step 6.2
    for each pdb with ending *cycle1.pdb:
    use src/utils/seq_frompdb.py
    '''
    mpnn_pdbs = glob.glob(f"{rf_dir}/*cycle*pdb")
    seq_dict = {'chainA': [], 
                'chainB': [],
                'coulen': [],
                'ljen': [],
                'mden': []}

    for it_f, pdb_it in enumerate(mpnn_pdbs):
        seq_new = seq_frompdb.get_seq_from_pdb(pdb_it)
        seq_dict['chainA'].append(seq_new[0])
        seq_dict['chainB'].append(seq_new[1])
        print(seq_new)
        seq_dict[it_f] = seq_new
        try:
            os.mkdir(f'{rf_dir}/chai_struct_{it_f}')
        except:
            pass
        chai_pred.fold_chai(
                seq_new[0],
                seq_new[1],
                f'{rf_dir}/chai_struct_{it_f}/temp.fasta',
                f'{rf_dir}/chai_struct_{it_f}',
                device=0
                )
        cif2pdb.cif2pdb(f'{rf_dir}/chai_struct_{it_f}/pred.model_idx_0.cif')

    '''
    step 7
    evaluate interaction energy for each generated structure
    '''
    for it_f, _ in enumerate(mpnn_pdbs):
        coul_en, lj_en = lie.lie(f'{rf_dir}/chai_struct_{it_f}/pred.model_idx_0.pdb',
                            f'{rf_dir}/chai_struct_{it_f}/fixed_0.pdb',
                            chaintarget,
                            rid_init,
                            rid_fin,
                            )
        md_en = md_energy.calculate_md_energy(
                            f'{rf_dir}/chai_struct_{it_f}/fixed_0.pdb',
                            )
        seq_dict['coulen'].append(coul_en)
        seq_dict['ljen'].append(lj_en)
        seq_dict['mden'].append(md_en)
    return seq_dict 




def run_pipeline(datafile,
                chai_dir,
                rf_dir,
                num_designs):

    data = pd.read_csv(datafile)
    device_ind = 0
    os.environ["CUDA_VISIBLE_DEVICES"] = str(device_ind)
    rf_script_path = "/lus/eagle/projects/datascience/avasan/RFDiffusionProject/RFdiffusion/scripts"
    se3_env = "/lus/eagle/projects/datascience/avasan/envs/SE3nv"
    dlbind_env = "/lus/eagle/projects/datascience/avasan/envs/proteinmpnn_binder_design"

    for rowit in range(1):#len(data)):
        # test with heavy_cdr3 without antigen
        chai_dir_it = f'{chai_dir}/{rowit}'
        rf_dir_it = f'{rf_dir}/{rowit}'
        try:
            os.mkdir(chai_dir_it)
        except:
            pass

        try:
            os.mkdir(rf_dir_it)
        except:
            pass

        cdrloop = "heavy_cdr3"
        chaintarget = "heavy_chain"

        seq_dict = run_pipeline_single(rowit,
                        data,
                        cdrloop,
                        chaintarget,
                        chai_dir_it,
                        rf_dir_it,
                        rf_script_path,
                        num_designs,
                        se3_env,
                        dlbind_env)

        df_seq = pd.DataFrame(seq_dict)
        df_seq.to_csv(f'seq_{rowit}.csv', index=False)




if __name__ =="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-i',
                        '--inputfil',
                        type=str,
                        help='input csv file with sequences')

    parser.add_argument('-C',
                        '--chaiout',
                        type=str,
                        help='directory to store chai output')

    parser.add_argument('-R',
                        '--rfout',
                        type=str,
                        help='directory to store rfdiff output')

    parser.add_argument('-N',
                        '--ndesigns',
                        type=str,
                        help='number of rfdiff designs')

    args = parser.parse_args()

    #'iedb_tables_pos_and_neg/lp_hl_np_p_ab_data.csv',

    try:
        os.mkdir(args.chaiout)
    except:
        pass

    try:
        os.mkdir(args.rfout)
    except:
        pass

    run_pipeline(args.inputfil,
                args.chaiout,
                args.rfout,
                args.ndesigns)
