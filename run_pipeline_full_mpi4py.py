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
                        cdrloop,
                        chaintarget,
                        chai_dir,
                        rf_dir,
                        residue_mapping,
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
    data_it_heavy = data_it['heavy_chain']
    data_it_light = data_it['light_chain']
    data_it_ant = data_it['antigen']

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

def chai_folding_seq(
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
    #os.environ["CUDA_VISIBLE_DEVICES"] = str(device_ind)
    os.environ["CHAI_DOWNLOADS_DIR"] = "/lus/eagle/projects/datascience/avasan/Software/CHAI1_downloads"

    rf_script_path = "/lus/eagle/projects/datascience/avasan/RFDiffusionProject/RFdiffusion/scripts"
    se3_env = "/lus/eagle/projects/datascience/avasan/envs/SE3nv"
    dlbind_env = "/lus/eagle/projects/datascience/avasan/envs/proteinmpnn_binder_design"

    for rowit in range(1):#len(data)):
        chai_dir_it = f'{chai_dir}/{rowit}'
        rf_dir_it = f'{rf_dir}/{rowit}'
        log_dir_it = f'{log_dir}/{rowit}'
        res_map_loc_it = f'{res_map_loc}/{rowit}'

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
