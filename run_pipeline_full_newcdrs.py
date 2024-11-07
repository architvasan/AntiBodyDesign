#import src.utils.id_cdrloop as cdrloop
import src.fold_ai.chai_pred as chai_pred
import os
import src.utils.cif2pdb as cif2pdb
import src.utils.fix_rfdiff_out as fix_rfdiff
import src.utils.seq_frompdb as seq_frompdb
import glob
import src.analyze.lie as lie
import src.rfdiffrun.partialdiff_loop as partial_loop
import src.rfdiffrun.fixrfpdb_seq as fixrfpdb_seq
import src.rfdiffrun.silenttools as silenttools
import src.rfdiffrun.dlbinder as dlbinder
import src.analyze.md_energy as md_energy
import src.utils.truncate as truncate
import src.utils.len_chains_pdb as len_chains
import pandas as pd
from mpi4py import MPI
import json
import numpy as np

'''
Single rank steps:
1. Load light + heavy chains into CHAI-1
2. Analyze md energy for generated structure
3. Truncate CHAI-1 structure to all residues within 15 A of CDR loop 

Parallel Steps:
0. Load in pandas dataframe
1. Determine sequencs for light + heavy chains + antigens
2. Determine resid numbers for CDR loop to target
4. Load truncated CHAI-1 structure into rfdiffusion 
    to generate new CDR loops with partial diff
4.1 Fix residues with known sequences.
4.2 Convert pdbs in directory to .silent file
5 Load rfdiffusion .silent file into dl binder design
6. Load dlbind design sequence into CHAI-1
6.1 Convert .silent files back to .pdb
8. Run a short MD/minimization simulation for each structure
9. Determine interaction energy + md energy from each simulation
'''

'''
Running in parallel using mpi4py
Just initialize MPI and 
set CUDA Visible device according to rank
'''

def initialize_mpi(gpu_per_node):
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()
    device = rank % gpu_per_node 
    return comm, size, rank, device


def chai_fold_pdb(
                    seq_a,
                    seq_b,
                    seq_c,
                    output_fasta,
                    output_dir,
                    ):
    chai_pred.fold_chai_body_ant(
            seq_a,
            seq_b,
            seq_c,
            output_fasta,
            output_dir,
            device=0,
            )
    cif2pdb.cif2pdb(f'{output_dir}/pred.model_idx_0.cif')
    return 

def energy_calcs(inp_pdb,
                 fixed_pdb,
                 chain_use,
                 rid_init,
                 rid_fin):

    ''' This following function 
        creates a fixed pdb 
        using pdbfixer
    '''
    coul_en, lj_en = lie.lie(inp_pdb,
                            fixed_pdb,
                            chain_use,
                            rid_init,
                            rid_fin,
                            )

    md_en = md_energy.calculate_md_energy(
                    fixed_pdb,
                    )
    return coul_en, lj_en, md_en

'''
Functions to truncate the structure using pymol2
'''
def truncate_functions(inp_pdb,
                       truncation_sel,
                       truncated_pdb,
                       resmap_outfile,
                       rid_init,
                       rid_fin):

    residue_mapping = truncate.truncate_pdb(
                                inp_pdb,
                                truncation_sel,
                                truncated_pdb,
                                )
    truncate.save_resmap(residue_mapping,
                                 resmap_outfile)

    rid_init_truncated = residue_mapping[('A', str(rid_init))]
    rid_fin_truncated = residue_mapping[('A', str(rid_fin))]

    return residue_mapping, rid_init_truncated, rid_fin_truncated

def rfdiffusion_generations(truncated_pdb,
                         rf_out,
                         rid_init,
                         rid_fin,
                         script_path,
                         se3_env,
                         device_ind,
                         num_designs,
                         ):
    len_A = len_chains.count_residues_in_chain(truncated_pdb,
                                                'A')
    len_B = len_chains.count_residues_in_chain(truncated_pdb,
                                                'B')
    len_C = len_chains.count_residues_in_chain(truncated_pdb,
                                                'C')
    '''
    Run RF Diffusion 
    Partial diffusion to diffuse new structure around the cdr loop
    (rid_init-rid_fin for chain A)
    need to load in separate conda env for se3_env here
    '''
    partial_loop.rfdiff_full(
              truncated_pdb,
              rf_out,
              rid_init,
              rid_fin,
              len_A,
              len_B,
              len_C,
              script_path,
              se3_env,
              device_ind,
              partial_steps=10,
              num_designs=num_designs,
              ) 
    '''
    step 4.15
    move any residue with > len(heavy_chain) to chain B + > len(heavy_chain+light_chain) to chain C
    '''
    for file in os.listdir(rf_out):
        print(file)
        if file.endswith(".pdb"):
            print(file)
            pdb_path = os.path.join(rf_out, file)
            fix_rfdiff.modify_chain_full(pdb_path,
                                         len_A,
                                         len_B,
                                        )
    return

def dlbinder_seq_gens(
                      local_pwd,
                      rf_out,
                      dlbind_env,
                      helperscripts,
                      silenttools_loc,
                      mpnn_fr,
                      device_ind,
                    ):
    
    '''
    fix residues in pdb that already have predefined sequence 
    i.e. not cdr loops 
    '''
    fixrfpdb_seq.fixpdb(rf_out,
                        dlbind_env,
                        helperscripts)
    
    '''convert pdb to silent file'''
    silenttools.pdb2silent(dlbind_env,
                           rf_out,
                           silenttools_loc)

    '''
    Run inverse folding with mpnn model
    '''
    dlbinder.protein_mpnn(dlbind_env,
                              mpnn_fr,
                              rf_out,
                              device_ind,
                        )

    '''
    convert silent file from mpnn back to pdb files
    '''
    silenttools.extractpdb(local_pwd,
                           dlbind_env,
                           rf_out,
                           silenttools_loc,
                           )

    return

def newcdr_folding(
                rf_out,
                rid_init,
                rid_fin,
                seq_A,
                seq_B,
                seq_C,
                seq_dict,
                ):
    '''
    step 6.2
    for each pdb with ending *cycle1.pdb:
    use src/utils/seq_frompdb.py
    '''
    mpnn_pdbs = glob.glob(f"{rf_out}/*cycle*pdb")

    for it_f, pdb_it in enumerate(mpnn_pdbs):
        '''
        get sequene for each pdb and plug into chai
        '''
        seq_new = seq_frompdb.get_seq_from_pdb(pdb_it)
        cdrnew_list = seq_new[0][rid_init:rid_fin]
        cdrnew = "".join(cdrnew_list)
        seq_A_new = seq_A[:rid_init] + cdrnew + seq_A[rid_fin:] 
        seq_dict['cdrseq'].append(cdrnew)
        try:
            os.mkdir(f'{rf_out}/chai_struct_{it_f}')
        except:
            pass

        chai_fold_pdb(
            seq_A_new,
            seq_B,
            seq_C,
            f'{rf_out}/chai_struct_{it_f}/temp.fasta',
            f'{rf_out}/chai_struct_{it_f}',
        )

        '''
        evaluate interaction + md energy for each generated structure
        '''
        coul_en, lj_en, md_en = energy_calcs(f'{rf_out}/chai_struct_{it_f}/pred.model_idx_0.pdb',
                                            f'{rf_out}/chai_struct_{it_f}/fixed_0.pdb',
                                            'A',
                                            rid_init,
                                            rid_fin)
        seq_dict['coulen'].append(coul_en)
        seq_dict['ljen'].append(lj_en)
        seq_dict['mden'].append(md_en)
        print(seq_dict)
    return seq_dict 


def full_pipeline_simple(df_db,
                         cdr_list,
                         chaintarget,
                         cdr_id,
                         output_dir,
                         rfdiff_scripts,
                         dlbind_helper_scripts,
                         silenttools_loc,
                         mpnn_fr,
                         local_pwd,
                         se3_env,
                         dl_bind_env,
                         designs_per_cdr=1,
                         device_ind = "0"):

    try:
        os.mkdir(output_dir)
    except:
        pass

    '''
    Set up dictionary we will write to for each of the generated cdr sequences
    '''
    dict_info = {'cdrseq': [],
                 'coulen': [],
                 'ljen': [],
                 'mden': [], 
                  }
    for rowit in range(1):
        '''
        1. Obtain info from original structure about:
        heavy chain
        light chain
        antigen
        cdr of interest
        rid_init_orig + rid_fin_orig for cdr
        '''
        try:
            os.mkdir(f'{output_dir}/{rowit}')
        except:
            pass

        df_it = df_db.loc[rowit]
        heavy_seq_orig = df_it['heavy_chain']
        light_seq_orig = df_it['light_chain']    
        ant_seq = df_it['antigen']
        orig_cdr_seq = df_it[cdr_id]
        
        if chaintarget == 'heavy_chain':
            rid_init_orig = heavy_seq_orig.find(orig_cdr_seq)
        elif chaintarget == 'light_chain':
            rid_init_orig = light_seq_orig.find(orig_cdr_seq)
        
        rid_fin_orig = rid_init_orig + len(orig_cdr_seq)
        

        
        for it_cdr, cdr in enumerate(cdr_list):
            '''
            2.0 create important directories in output_dir
            '''
            try:
                os.mkdir(f'{output_dir}/{rowit}/{it_cdr}')
                os.mkdir(f'{output_dir}/{rowit}/{it_cdr}/chaiout')
                os.mkdir(f'{output_dir}/{rowit}/{it_cdr}/rfout')
                os.mkdir(f'{output_dir}/{rowit}/{it_cdr}/resmaps')
            except:
                pass

            '''
            2. graft generated cdr into new sequence and generate structure
            set seq_new (with new cdr pasted in) whether heavy or light chain is target
            fold sequences using chai (for both heavy and light chains)
            Keep sequence we want to affect with rf diffusion as the first sequence
            Write to ...chaiout/pred.model_idx_0.pdb
            '''
            if chaintarget == 'heavy_chain':
                seq_new = heavy_seq_orig[:rid_init_orig] + cdr + heavy_seq_orig[rid_fin_orig:]
                chai_fold_pdb(
                    seq_new,
                    light_seq_orig,
                    ant_seq,
                    f'{output_dir}/{rowit}/{it_cdr}/chaiout/temp.fasta',
                    f'{output_dir}/{rowit}/{it_cdr}/chaiout',
                    device='0',
                    )
            elif chaintarget == 'light_chain':
                seq_new = light_seq_orig[:rid_init_orig] + cdr + light_seq_orig[rid_fin_orig:]
                chai_fold_pdb(
                    seq_new,
                    heavy_seq_orig,
                    ant_seq,
                    f'{output_dir}/{rowit}/{it_cdr}/chaiout/temp.fasta',
                    f'{output_dir}/{rowit}/{it_cdr}/chaiout',
                    device='0',
                    )

            '''
            3. evaluate coul + lj + md energy for new structure
            '''
            rid_init_new = rid_init_orig
            rid_fin_new = rid_init_new + len(cdr)
            chai_gen_pdb = f'{output_dir}/{rowit}/{it_cdr}/chaiout/pred.model_idx_0.pdb'
            chai_fixed_pdb = f'{output_dir}/{rowit}/{it_cdr}/chaiout/fixed.pdb'

            coul_en, lj_en, md_en = energy_calcs(chai_gen_pdb,
                                                chai_fixed_pdb,
                                                'A',
                                                rid_init_new,
                                                rid_fin_new)

            '''
            write cdr seq + energies to dict
            '''
            dict_info['cdrseq'].append(cdr)
            dict_info['coulen'].append(coul_en)
            dict_info['ljen'].append(lj_en)
            dict_info['mden'].append(md_en)
            if it_cdr % 10 ==0:
                print(dict_info)
            
            '''
            RF Diffusion stuff
            '''
            '''
            4. Truncate chaigen structure
                to within 15 A of selection
            '''
            chai_truncated_pdb = f'{output_dir}/{rowit}/{it_cdr}/chaiout/pred_0_truncated.pdb'
            truncation_sel = f'chain A and resi {rid_init_new}-{rid_fin_new}',
            resmap_outfile = f'{output_dir}/{rowit}/{it_cdr}/resmaps/resmap_0.pkl'

            residue_mapping,\
                  rid_init_truncated,\
                      rid_fin_truncated = truncate_functions(
                                            chai_gen_pdb,
                                            truncation_sel,
                                            chai_truncated_pdb, 
                                            resmap_outfile,
                                            rid_init_new,
                                            rid_fin_new,
                                            )

            '''
            5. Run rfdiffusion to generate n designs in rfout location
               Partial difusion around cdr loop
            '''
            rfdiffusion_generations(
                         chai_truncated_pdb,
                         f'{output_dir}/{rowit}/{it_cdr}/rfout',
                         rid_init_truncated,
                         rid_fin_truncated,
                         rfdiff_scripts,
                         se3_env,
                         device_ind,
                         designs_per_cdr,
                         )

            '''
            6. Run dlbinder mpnn to obtain the inverse folded structure 
               convert glycines in cdr region to realistic amino acids
            '''
            dlbinder_seq_gens(
                      local_pwd,
                      f'{output_dir}/{rowit}/{it_cdr}/rfout',
                      dl_bind_env,
                      dlbind_helper_scripts,
                      silenttools_loc,
                      mpnn_fr,
                      device_ind,
                    )
            '''
            7.step 6.2
              for each pdb generated from mpnn with ending *cycle1.pdb:
              use src/utils/seq_frompdb.py to get amino acid sequence for cdrloop 
            '''
            if chaintarget=='heavy_chain':
                dict_info = newcdr_folding(
                                        f'{output_dir}/{rowit}/{it_cdr}/rfout',
                                        rid_init_truncated,
                                        rid_fin_truncated,
                                        seq_new,
                                        light_seq_orig,
                                        ant_seq,
                                        dict_info,
                                        )
            elif chaintarget=='light_chain':
                dict_info = newcdr_folding(
                                        f'{output_dir}/{rowit}/{it_cdr}/rfout',
                                        rid_init_truncated,
                                        rid_fin_truncated,
                                        seq_new,
                                        heavy_seq_orig,
                                        ant_seq,
                                        dict_info,
                                        )
    return dict_info 

def run_pipeline(db_file,
                 cdrlistfile,
                 mpi_yes,
                 chaintarget,
                 cdr_id,
                 output_dir,
                 designs_per_cdr,
                 local_pwd,
                 rfdiff_scripts,
                 dlbind_helper_scripts,
                 silenttools_loc,
                 mpnn_fr, 
                 se3_env,
                 dl_bind_env,
                 gpu_per_node=4):

    df_db = pd.read_csv(db_file)
    with open(cdrlistfile, 'r') as file:
        cdrlist = json.load(file) 
    try:
        os.mkdir(output_dir)
    except:
        pass
    if mpi_yes:
        comm, size, rank, device = initialize_mpi(gpu_per_node)
        cdrlist_rank = np.array_split(cdrlist, int(size))[int(rank)]
        try:
            os.mkdir(f'{output_dir}/{rank}')
        except:
            pass

        os.environ["CUDA_VISIBLE_DEVICES"] = str(device)
        os.environ["CHAI_DOWNLOADS_DIR"] = "/lus/eagle/projects/datascience/avasan/Software/CHAI1_downloads"

        dict_info = full_pipeline_simple(df_db,
                         cdrlist_rank,
                         chaintarget,
                         cdr_id,
                         f'{output_dir}/{rank}',
                         local_pwd,
                         rfdiff_scripts,
                         dlbind_helper_scripts,
                         silenttools_loc,
                         mpnn_fr,
                         se3_env,
                         dl_bind_env,
                         designs_per_cdr=designs_per_cdr,
                         device_ind = device)

        df_info = pd.DataFrame(dict_info)
        df_info.to_csv(f'{output_dir}/{rank}/df_cdrs_energies.csv')
    else:
        device = 0
        dict_info = full_pipeline_simple(df_db,
                         cdrlist_rank,
                         chaintarget,
                         cdr_id,
                         f'{output_dir}',
                         local_pwd,
                         rfdiff_scripts,
                         dlbind_helper_scripts,
                         silenttools_loc,
                         mpnn_fr,
                         se3_env,
                         dl_bind_env,
                         designs_per_cdr=1,
                         device_ind = device)
        df_info = pd.DataFrame(dict_info)
        df_info.to_csv(f'{output_dir}/df_cdrs_energies.csv')
    return dict_info

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-df',
                        '--dbfil',
                        type=str,
                        help='input csv file with sequences')

    parser.add_argument('-cf',
                        '--cdrfil',
                        type=str,
                        help='file with cdr seqs')
    
    parser.add_argument('-C',
                        '--chaintarget',
                        type=str,
                        help='file with cdr seqs')

    parser.add_argument('-CD',
                        '--cdrid',
                        type=str,
                        help='file with cdr seqs')

    parser.add_argument('-N',
                        '--ndesigns',
                        type=int,
                        required=False,
                        default=1,
                        help='number of rfdiff designs per cdr')   
    
    parser.add_argument('-O',
                        '--outdir',
                        type=str,
                        required=True,
                        help='output directory')   

    parser.add_argument('-m',
                        '--mpiyes',
                        type=bool,
                        help='use mpi or not?')

    parser.add_argument('-G',                              
                         '--gpunum',
                         type=int,
                         required=False,
                         default=4,
                         help='number of gpus on a single node')

    args = parser.parse_args()

    rfdiff_scripts = ""
    rfdiff_scripts = "/lus/eagle/projects/datascience/avasan/RFDiffusionProject/RFdiffusion/scripts"
    se3_env = "/lus/eagle/projects/datascience/avasan/envs/SE3nv"
    dl_bind_env = "/lus/eagle/projects/datascience/avasan/envs/proteinmpnn_binder_design"
    dlbind_helper_scripts = '/eagle/projects/datascience/avasan/RFDiffusionProject/dl_binder_design/helper_scripts'
    silenttools_loc = '/eagle/projects/datascience/avasan/RFDiffusionProject/dl_binder_design/include/silent_tools'
    mpnn_fr = '/eagle/projects/datascience/avasan/RFDiffusionProject/dl_binder_design/mpnn_fr'
    local_pwd = '/eagle/datascience/avasan/Simulations/AntiBodyDesign'
    
    run_pipeline(args.dbfil,
                 args.cdrfil,
                 args.mpiyes,
                 args.chaintarget,
                 args.cdrid,
                 args.outdir,
                 args.ndesigns,
                 local_pwd,
                 rfdiff_scripts,
                 dlbind_helper_scripts,
                 silenttools_loc,
                 mpnn_fr, 
                 se3_env,
                 dl_bind_env,
                 gpu_per_node=args.gpunum)


if False:
    '''
    Single rank steps:
    1. Load light + heavy chains into CHAI-1
    2. Analyze md energy for generated structure
    3. Truncate CHAI-1 structure to all residues within 15 A of CDR loop 

    Parallel Steps:
    0. Load in pandas dataframe
    1. Determine sequencs for light + heavy chains + antigens
    2. Determine resid numbers for CDR loop to target
    4. Load truncated CHAI-1 structure into rfdiffusion 
        to generate new CDR loops with partial diff
    4.1 Fix residues with known sequences.
    4.2 Convert pdbs in directory to .silent file
    5 Load rfdiffusion .silent file into dl binder design
    6. Load dlbind design sequence into CHAI-1
    6.1 Convert .silent files back to .pdb
    8. Run a short MD/minimization simulation for each structure
    9. Determine interaction energy + md energy from each simulation
    '''

    def run_pipeline_single(
                            rowit,
                            data,
                            cdrlist,
                            cdrloop,
                            chaintarget,
                            chai_dir,
                            rf_dir,
                            log_dir,
                            residue_mapping_loc,
                            rf_script_path,
                            num_designs,
                            se3_env,
                            dlbind_env,
                            device_ind
                            ):

        seq_dict = {'it': [],
                    'chainA': [], 
                    'chainB': [],
                    'antigen': [],
                    'cdrseq': [],
                    'coulen': [],
                    'ljen': [],
                    'mden': [],
                    'struct': []}

        '''
        step 1
        '''
        data_it = data.loc[rowit]
        heavy_seq = data_it['heavy_chain']
        light_seq = data_it['light_chain']
        ant_seq = data_it['antigen']
        cdrloop_seq = cdrlist[rowit]

        '''
        step 2
        '''
        data_it_target = data_it[chaintarget]
        data_it_patt = data_it[cdrloop]

        rid_init = data_it_target.find(data_it_patt)
        rid_fin = rid_init + len(data_it_patt)

        if chaintarget == 'heavy_chain':
            chain_use = 'A'
        elif chaintarget == 'light_chain':
            chain_use = 'B'

        seq_dict_init, truncated_res_map = chai_folding_seq_new_cdr(
                                                            heavy_seq,
                                                            light_seq,
                                                            ant_seq,
                                                            cdrloop_seq,
                                                            chaintarget,
                                                            rid_init,
                                                            rid_fin,
                                                            chai_dir,
                                                            map_dir,
                                                            rank)
        seq_init_df = pd.DataFrame([seq_dict_init])
        seq_init_df.to_csv(f'{log_dir_it}/seq_initial_{rowit}_{cdrloop}.csv', index=False)




        '''
        step 3
        Construct and execute the command
        Use truncated system here
        Use preloaded residue mapping dictionary
        * Assume for now that we are dealing with heavy chain.
        Will modify for light chain later
        '''

        print("residue mapping is:")
        print(residue_mapping)

        rid_init_truncated = residue_mapping[('A', str(rid_init))]
        rid_fin_truncated = residue_mapping[('A', str(rid_fin))]

        len_heavy = len_chains.count_residues_in_chain(f'{chai_dir}/pred_0_truncated.pdb', 'A')
        len_light = len_chains.count_residues_in_chain(f'{chai_dir}/pred_0_truncated.pdb', 'B')
        len_ant = len_chains.count_residues_in_chain(f'{chai_dir}/pred_0_truncated.pdb', 'C')

        print(chai_dir)
        if True:
            partial_loop.rfdiff_full(
                  f'{chai_dir}/pred_0_truncated.pdb',
                  rf_dir,
                  rid_init_truncated,
                  rid_fin_truncated,
                  len_heavy,
                  len_light,
                  len_ant,
                  rf_script_path,
                  se3_env,
                  device_ind,
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
                                      rf_dir,
                                      device_ind,
                                )

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

        for it_f, pdb_it in enumerate(mpnn_pdbs):
            seq_new = seq_frompdb.get_seq_from_pdb(pdb_it)
            cdrnew_list = seq_new[0][rid_init_truncated:rid_fin_truncated]
            cdrnew = "".join(cdrnew_list)
            seq_A_new = data_it_heavy[:rid_init] + cdrnew + data_it_heavy[rid_fin:] 
            seq_dict['it'].append(it_f + 1)
            seq_dict['chainA'].append(seq_A_new)
            seq_dict['chainB'].append(data_it_light)
            seq_dict['antigen'].append(data_it_ant)
            seq_dict['cdrseq'].append(cdrnew)

            try:
                os.mkdir(f'{rf_dir}/chai_struct_{it_f}')
            except:
                pass

            chai_pred.fold_chai_body_ant(
                                seq_A_new,
                                data_it_light,
                                data_it_ant,
                                f'{rf_dir}/chai_struct_{it_f}/temp.fasta',
                                f'{rf_dir}/chai_struct_{it_f}',
                                device=0
                                )

            cif2pdb.cif2pdb(f'{rf_dir}/chai_struct_{it_f}/pred.model_idx_0.cif')
            seq_dict['struct'].append(f'{rf_dir}/chai_struct_{it_f}/pred.model_idx_0.pdb')

            '''
            step 7
            evaluate interaction + md energy for each generated structure
            '''
            coul_en, lj_en = lie.lie(f'{rf_dir}/chai_struct_{it_f}/pred.model_idx_0.pdb',
                                f'{rf_dir}/chai_struct_{it_f}/fixed_0.pdb',
                                chain_use,
                                rid_init,
                                rid_fin,
                                )

            md_en = md_energy.calculate_md_energy(
                                f'{rf_dir}/chai_struct_{it_f}/fixed_0.pdb',
                                )

            seq_dict['coulen'].append(coul_en)
            seq_dict['ljen'].append(lj_en)
            seq_dict['mden'].append(md_en)
            print(seq_dict)
        return seq_dict 

    def chai_folding_seq_new_cdr(
                        heavy_seq,
                        light_seq,
                        ant_seq,
                        cdrloop_seq,
                        chaintarget,
                        rid_init,
                        rid_fin,
                        chai_dir,
                        map_dir,
                        rank):
                        #comm):

        heavy_seq_new = heavy_seq[:rid_init] + cdrloop_seq + heavy_seq[rid_fin:] 
        chai_pred.fold_chai_body_ant(
                heavy_seq_new,
                light_seq,
                ant_seq,
                f'{chai_dir}/temp.fasta',
                chai_dir,
                device=0
                )
        '''
        use cif2pdb util to convert cif to pdb
        '''
        cif2pdb.cif2pdb(f'{chai_dir}/pred.model_idx_0.cif')

        if chaintarget == 'heavy_chain':
             chain_use = 'A'
        elif chaintarget == 'light_chain':
             chain_use = 'B'

        rid_init_new = heavy_seq_new.find(cdrloop_seq)
        rid_fin_new = rid_init + len(cdrloop_seq)

        coul_en, lj_en = lie.lie(f'{chai_dir}/pred.model_idx_0.pdb',
                              f'{chai_dir}/fixed_0.pdb',
                              chain_use,
                              rid_init_new,
                              rid_fin_new,
                              )

        md_en = md_energy.calculate_md_energy(
                    f'{chai_dir}/fixed_0.pdb',
                    )

        residue_mapping = truncate.truncate_pdb(
                                f'{chai_dir}/pred.model_idx_0.pdb',
                                f'chain A and resi {rid_init_new}-{rid_fin_new}',
                                f'{chai_dir}/pred_0_truncated.pdb',
                                )
        truncate.save_resmap(residue_mapping, f'{map_dir}/resmap_0.pkl')

        seq_dict = {'it': rowit,
             'chainA': data_it_heavy, 
             'chainB': data_it_light,
             'antigen': data_it_ant,
             'cdrseq': data_it_patt,
             'coulen': coul_en,
             'ljen': lj_en,
             'mden': md_en,
             'struct': f'{chai_dir}/pred.model_idx_0.pdb'}         
            #comm.Barrier()
        return seq_dict, residue_mapping





    def _chai_folding_seq(
                        data,
                        rowit,
                        chai_dir,
                        chaintarget,
                        cdrloop,
                        map_dir,
                        rank):
                        #comm):

        if True:
            data_it = data.loc[rowit]
            data_it_heavy = data_it['heavy_chain']
            data_it_light = data_it['light_chain']
            data_it_ant = data_it['antigen']
            chai_pred.fold_chai_body_ant(
                                data_it_heavy,
                                data_it_light,
                                data_it_ant,
                                f'{chai_dir}/temp.fasta',
                                chai_dir,
                                device=0
                                )
            '''
            use cif2pdb util to convert cif to pdb
            '''
            cif2pdb.cif2pdb(f'{chai_dir}/pred.model_idx_0.cif')

            if chaintarget == 'heavy_chain':
                 chain_use = 'A'
            elif chaintarget == 'light_chain':
                 chain_use = 'B'
            data_it_target = data_it[chaintarget]
            data_it_patt = data_it[cdrloop]

            rid_init = data_it_target.find(data_it_patt)
            rid_fin = rid_init + len(data_it_patt)

            coul_en, lj_en = lie.lie(f'{chai_dir}/pred.model_idx_0.pdb',
                                  f'{chai_dir}/fixed_0.pdb',
                                  chain_use,
                                  rid_init,
                                  rid_fin,
                                  )

            md_en = md_energy.calculate_md_energy(
                                f'{chai_dir}/fixed_0.pdb',
                                )

            residue_mapping = truncate.truncate_pdb(
                                    f'{chai_dir}/pred.model_idx_0.pdb',
                                    f'chain A and resi {rid_init}-{rid_fin}',
                                    f'{chai_dir}/pred_0_truncated.pdb',
                                    )
            truncate.save_resmap(residue_mapping, f'{map_dir}/resmap_0.pkl')

            seq_dict = {'it': rowit,
                 'chainA': data_it_heavy, 
                 'chainB': data_it_light,
                 'antigen': data_it_ant,
                 'cdrseq': data_it_patt,
                 'coulen': coul_en,
                 'ljen': lj_en,
                 'mden': md_en,
                 'struct': f'{chai_dir}/pred.model_idx_0.pdb'}         
            #comm.Barrier()
        return seq_dict, residue_mapping

    def run_pipeline_full(
                    datafile,
                    cdrlistfile,
                    chai_dir,
                    rf_dir,
                    log_dir,
                    res_map_loc,
                    num_designs,
                    gpu_per_node=4,
                    fold_init = False):

        '''
        Initialize mpi rank + device index
        '''

        data = pd.read_csv(datafile)
        with open(cdrlistfile, 'r') as file:
            cdrlist = json.load(file)

        #os.environ["CUDA_VISIBLE_DEVICES"] = str(device_ind)
        os.environ["CHAI_DOWNLOADS_DIR"] = "/lus/eagle/projects/datascience/avasan/Software/CHAI1_downloads"

        rf_script_path = "/lus/eagle/projects/datascience/avasan/RFDiffusionProject/RFdiffusion/scripts"
        se3_env = "/lus/eagle/projects/datascience/avasan/envs/SE3nv"
        dlbind_env = "/lus/eagle/projects/datascience/avasan/envs/proteinmpnn_binder_design"

        for rowit in range(1):#len(data)):
            for cdrit in range(len(cdrlist)):
                chai_dir_it = f'{chai_dir}/{rowit}_{cdrit}'
                rf_dir_it = f'{rf_dir}/{rowit}_{cdrit}'
                log_dir_it = f'{log_dir}/{rowit}'
                res_map_loc_it = f'{res_map_loc}/{rowit}_{cdrit}'

                try:
                    os.mkdir(res_map_loc_it)
                except:
                    pass

                try:
                    os.mkdir(chai_dir_it)
                except:
                    pass

                try:
                    os.mkdir(rf_dir_it)
                except:
                    pass

                try:
                    os.mkdir(log_dir_it)
                except:
                    pass

                cdrloop = "heavy_cdr3"
                chaintarget = "heavy_chain"

                print(fold_init)
                print(chai_dir_it)
                if fold_init == True:
                    rank = 0

                    seq_dict_init, truncated_res_map = chai_folding_seq(
                                            data,
                                            rowit,
                                            chai_dir_it,
                                            chaintarget,
                                            cdrloop,
                                            res_map_loc_it,
                                            rank)
                                            #comm)
                    seq_init_df = pd.DataFrame([seq_dict_init])
                    seq_init_df.to_csv(f'{log_dir_it}/seq_initial_{rowit}_{cdrloop}.csv', index=False)

                else:
                    truncated_res_map = truncate.open_resmap(f"{res_map_loc_it}/resmap_{rowit}.pkl")


                seq_dict = run_pipeline_single(
                                rowit,
                                data,
                                cdrloop,
                                chaintarget,
                                chai_dir_it,
                                rf_dir_it,
                                truncated_res_map,
                                rf_script_path,
                                num_designs,
                                se3_env,
                                dlbind_env,
                                device_ind=0)

                print(seq_dict)
                df_seq = pd.DataFrame(seq_dict)
                df_seq.to_csv(f'{log_dir_it}/seq_{rowit}_{cdrloop}.csv', index=False)

    def run_pipeline_parallel(
                    datafile,
                    chai_dir,
                    rf_dir,
                    log_dir,
                    res_map_loc,
                    num_designs,
                    gpu_per_node=4,
                    fold_init = False):

        '''
        Initialize mpi rank + device index
        '''

        comm, size, rank, device_ind =  initialize_mpi(gpu_per_node)
        print(f"rank {rank} of {size}")
        data = pd.read_csv(datafile)
        print(device_ind)
        #os.environ["CUDA_VISIBLE_DEVICES"] = str(device_ind)
        os.environ["CHAI_DOWNLOADS_DIR"] = "/lus/eagle/projects/datascience/avasan/Software/CHAI1_downloads"

        rf_script_path = "/lus/eagle/projects/datascience/avasan/RFDiffusionProject/RFdiffusion/scripts"
        se3_env = "/lus/eagle/projects/datascience/avasan/envs/SE3nv"
        dlbind_env = "/lus/eagle/projects/datascience/avasan/envs/proteinmpnn_binder_design"

        for rowit in range(1):#len(data)):
            chai_dir_it = f'{chai_dir}/{rowit}'
            rf_dir_it = f'{rf_dir}/{rank}_{rowit}'
            log_dir_it = f'{log_dir}/{rank}_{rowit}'

            try:
                os.mkdir(chai_dir_it)
            except:
                pass

            try:
                os.mkdir(rf_dir_it)
            except:
                pass

            try:
                os.mkdir(log_dir_it)
            except:
                pass

            cdrloop = "heavy_cdr3"
            chaintarget = "heavy_chain"

            print(fold_init)
            if fold_init == True and rank == 0:
                seq_dict_init, truncated_res_map = chai_folding_seq(
                                        data,
                                        rowit,
                                        chai_dir_it,
                                        chaintarget,
                                        cdrloop,
                                        res_map_loc,
                                        rank,
                                        comm)
                seq_init_df = pd.DataFrame([seq_dict_init])
                seq_init_df.to_csv(f'{log_dir_it}/seq_initial_{rowit}_{cdrloop}.csv', index=False)

            else:
                truncated_res_map = truncate.open_resmap(f"{res_map_loc}/resmap_{rowit}.pkl")


            seq_dict = run_pipeline_single(
                            rowit,
                            data,
                            cdrloop,
                            chaintarget,
                            chai_dir_it,
                            rf_dir_it,
                            truncated_res_map,
                            rf_script_path,
                            num_designs,
                            se3_env,
                            dlbind_env,
                            device_ind=0)

            print(seq_dict)
            df_seq = pd.DataFrame(seq_dict)
            df_seq.to_csv(f'{log_dir_it}/seq_{rowit}_{cdrloop}.csv', index=False)

    if __name__ == "__main__":
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

        parser.add_argument('-L',
                            '--logout',
                            type=str,
                            help='directory to store log info (seqs, energies)')

        parser.add_argument('-M',
                            '--mapdir',
                            type=str,
                            help='directory where resmaps from truncation are')

        parser.add_argument('-N',
                            '--ndesigns',
                            type=int,
                            required=False,
                            default=10,
                            help='number of rfdiff designs')

        parser.add_argument('-G',                              
                             '--gpunum',
                             type=int,
                             required=False,
                             default=4,
                             help='number of gpus on a single node')

        parser.add_argument('-F',                              
                             '--foldinit',
                             type=bool,
                             required=False,
                             default=False,
                             help='should we fold the initial sequence (True) or is it prefolded? (False)')

        args = parser.parse_args()

        try:
            os.mkdir(args.chaiout)
        except:
            pass

        try:
            os.mkdir(args.rfout)
        except:
            pass

        try:                     
            os.mkdir(args.logout)
        except:
            pass

        try:
            os.mkdir(args.mapdir)
        except:
            pass
        run_pipeline_full(
                        args.inputfil,
                        args.chaiout,
                        args.rfout,
                        args.logout,
                        args.mapdir,
                        args.ndesigns,
                        gpu_per_node=args.gpunum,
                        fold_init = True)#args.foldinit)

    #'/eagle/projects/datascience/avasan/RFDiffusionProject/dl_binder_design/helper_scripts',
    #'/eagle/projects/datascience/avasan/RFDiffusionProject/dl_binder_design/include/silent_tools',
    #'/eagle/projects/datascience/avasan/RFDiffusionProject/dl_binder_design/mpnn_fr'