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
6. Load dlbind design into af2_predict
    get: new sequence + pae interaction score
6.1 Convert .silent files back to .pdb
7. Determine pae interaction score to filter out poorly designed structures
8. Run a short MD simulation for each structure
9. Determine average interaction energy from each simulation
'''

def run_pipeline_single(rowit,
                        data,
                        cdrloop,
                        chaintarget,
                        chai_dir,
                        rf_dir,
                        rf_script_path,
                        se3_env,
                        dlbind_env):
    '''
    step 1
    '''
    data_it = data.loc[rowit]
    data_it_heavy = data_it['heavy_chain']
    data_it_light = data_it['light_chain']

    '''
    step 2
    '''
    #data_it_ant = data_it['antigen']
    data_it_target = data_it[chaintarget]
    data_it_patt = data_it[cdrloop]
    rid_init = data_it_target.find(data_it_patt)
    rid_fin = rid_init + len(data_it_patt)

    if False:
        ### testing if cdrloop is correct
        data_it_cdr_test = list(data_it_target)[rid_init:rid_fin]
        print(''.join(data_it_cdr_test))
        print(data_it_patt)

    if False:
        '''
        step 3
        '''
        chai_pred.fold_chai(
                            data_it_heavy,
                            data_it_light,
                            f'{chai_dir}/temp.fasta',
                            chai_dir,
                            device=0
                            )

        ### use cif2pdb util to convert cif to pdb
        cif2pdb.cif2pdb(f'{chai_dir}/pred.model_idx_0.cif')

    '''
    step 4
    Construct and execute the command
    * Assume for now that we are dealing with heavy chain.
    Will modify for light chain later
    '''

    len_heavy = len(list(data_it_heavy))
    len_light = len(list(data_it_light))

    '''
    command = f"""module use /soft/modulefiles &&\
                module load conda &&\
                conda activate {se3_env}\
                && {rf_script_path}/run_inference.py\
                inference.output_prefix={rf_dir}/pred\
                inference.input_pdb={chai_dir}/pred.model_idx_0.pdb\
                'contigmap.contigs=["{len_heavy}-{len_heavy}/0 {len_light}-{len_light}"]'\
                diffuser.partial_T=10\
                'contigmap.provide_seq=[0-{rid_init-1},{rid_fin}-{len_heavy},{len_heavy+1}-{len_heavy+len_light}]'\
                inference.num_designs=1\
                """
    '''

    command = f"""module use /soft/modulefiles &&\
                module load conda &&\
                conda activate {se3_env}\
                && {rf_script_path}/run_inference.py\
                inference.output_prefix={rf_dir}/pred\
                inference.input_pdb={chai_dir}/pred.model_idx_0.pdb\
                'contigmap.contigs=[A1-{rid_init-1}/{rid_fin - rid_init}-{rid_fin - rid_init}/A{rid_fin}-{len_heavy}/0 B1-{len_light}]'\
                diffuser.partial_T=10\
                'contigmap.provide_seq=[0-{rid_init-1},{rid_fin}-{len_heavy},{len_heavy+1}-{len_heavy+len_light}]'\
                inference.num_designs=1\
                """


    if True:
        os.system(command)

    '''
    step 4.15
    move any residue with > len(heavy_chain) to chain B
    '''
    for file in os.listdir(rf_dir):
        print(file)
        if file.endswith(".pdb"):
            print(file)
            pdb_path = os.path.join(rf_dir, file)
            fix_rfdiff.modify_chain(pdb_path,
                                    len_heavy,
                                    )

    if True:
        '''
        step 4.1
        fix residues in pdb 
        '''
        command = f"""\
                  module use /soft/modulefiles &&\
                  module load conda &&\
                  conda activate /lus/eagle/projects/datascience/avasan/envs/proteinmpnn_binder_design &&\
                  python /eagle/projects/datascience/avasan/RFDiffusionProject/dl_binder_design/helper_scripts/addFIXEDlabels.py\
                  --pdbdir  {rf_dir}\
                  --trbdir {rf_dir}\
                  --verbose\
                  """
        print(command)
        os.system(command)

        '''
        step 4.2
        convert pdb to silent
        '''
        command = f"""\
                  module use /soft/modulefiles &&\
                  module load conda &&\
                  conda activate /lus/eagle/projects/datascience/avasan/envs/proteinmpnn_binder_design &&\
                  /eagle/projects/datascience/avasan/RFDiffusionProject/dl_binder_design/include/silent_tools/silentfrompdbs\
                  {rf_dir}/*.pdb >\
                  {rf_dir}/rfout.silent\
                  """

        print(command)
        os.system(command)

        '''
        step 5
        '''
        command = f"""\
                  module use /soft/modulefiles &&\
                  module load conda &&\
                  conda activate /lus/eagle/projects/datascience/avasan/envs/proteinmpnn_binder_design &&\
                  /eagle/projects/datascience/avasan/RFDiffusionProject/dl_binder_design/mpnn_fr/dl_interface_design.py\
                  -silent {rf_dir}/rfout.silent\
                  -outsilent {rf_dir}/mpnnout.silent\
                  -checkpoint_name {rf_dir}/checkpoint_mpnn.dat\
                 """
        print(command)
        os.system(command)

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
        command = f"""\
                  module use /soft/modulefiles &&\
                  module load conda &&\
                  cd {rf_dir} &&\
                  conda activate /lus/eagle/projects/datascience/avasan/envs/proteinmpnn_binder_design &&\
                  /eagle/projects/datascience/avasan/RFDiffusionProject/dl_binder_design/include/silent_tools/silentextract\
                  mpnnout.silent &&\
                  cd /eagle/datascience/avasan/Simulations/Antibody_Design\
                  """
        print(command)
        os.system(command)

        '''
        step 6.2
        for each pdb with ending *cycle1.pdb:
        use src/utils/seq_frompdb.py
        '''
        mpnn_pdbs = glob.glob(f"{rf_dir}/*cycle*pdb")
        seq_dict = {'chainA':[], 'chainB':[], 'coulen':[], 'ljen':[]}
        for it_f, pdb_it in enumerate(mpnn_pdbs):
            seq_new = get_seq_from_pdb(pdb_it)
            seq_dict['chainA'].append(seq_new[0])
            seq_dict['chainB'].append(seq_new[1])
            print(seq_new)
            seq_dict[it_f] = seq_new
            os.mkdir(f'{rf_dir}/chai_struct_{it_f}')
            chai_pred.fold_chai(
                    seq_new[0],
                    seq_new[1],
                    f'{rf_dir}/chai_struct_{it_f}/temp.fasta',
                    f'{rf_dir}/chai_struct_{it_f}',
                    device=0
                    )

        '''
        step 7
        evaluate interaction energy for each generated structure
        '''
        for it_f, _ in enumerate(mpnn_pdbs):
            coul_en, lj_en = lie(f'{rf_dir}/chai_struct_{it_f}/pred.model_idx_0.pdb',
                                f'{rf_dir}/chai_struct_{it_f}/fixed_0.pdb',
                                chaintarget,
                                rid_init,
                                rid_fin,
                                )
            seq_dict['coulen'].append(coul_en)
            seq_dict['lj_en'].append(lj_en)

        if False:
            '''
            step 6
            '''
            command = f"""\
                      module use /soft/modulefiles &&\
                      module load conda &&\
                      conda activate /lus/eagle/projects/datascience/avasan/envs/af2_binder_design &&\
                      /eagle/projects/datascience/avasan/RFDiffusionProject/dl_binder_design/af2_initial_guess/predict.py\
                      -silent {rf_dir}/mpnnout.silent\
                      -outsilent {rf_dir}/af2out.silent\
                      -checkpoint_name {rf_dir}/checkpoint_af2.dat\
                     """
            print(command)
            os.system(command)
            
            '''
            step 6.1
            '''
            command = f"""\
                      module use /soft/modulefiles &&\
                      module load conda &&\
                      cd {rf_dir} &&\
                      conda activate /lus/eagle/projects/datascience/avasan/envs/proteinmpnn_binder_design &&\
                      /eagle/projects/datascience/avasan/RFDiffusionProject/dl_binder_design/include/silent_tools/silentextract\
                      af2out.silent &&\
                      cd /eagle/datascience/avasan/Simulations/Antibody_Design\
                      """

            print(command)
            os.system(command)
            
            '''
            step 6.2
            '''
            return seq_new, coul_en, lj_en


def run_pipeline(datafile,
                chai_dir,
                rf_dir):
    data = pd.read_csv(datafile)
    device_ind = 0
    os.environ["CUDA_VISIBLE_DEVICES"] = str(device_ind)
    rf_script_path = "/lus/eagle/projects/datascience/avasan/RFDiffusionProject/RFdiffusion/scripts"
    se3_env = "/lus/eagle/projects/datascience/avasan/envs/SE3nv"
    dlbind_env = "/lus/eagle/projects/datascience/avasan/envs/af2_binder_design"

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

        run_pipeline_single(rowit,
                        data,
                        cdrloop,
                        chaintarget,
                        chai_dir_it,
                        rf_dir_it,
                        rf_script_path,
                        se3_env,
                        dlbind_env)

try:
    os.mkdir('chai_output')
except:
    pass
try:
    os.mkdir('rf_output')
except:
    pass

run_pipeline('iedb_tables_pos_and_neg/lp_hl_np_p_ab_data.csv',
            'chai_output',
            'rf_output')
