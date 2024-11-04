import os

def rfdiff(input_pdb,
           rf_dir,
           rid_init,
           rid_fin,
           len_heavy,
           len_light,
           rf_script_path,
           se3_env,
           partial_steps=10,
           num_designs=1,
           ):

    command = f"""module use /soft/modulefiles &&\
            module load conda &&\
            conda activate {se3_env}\
            && {rf_script_path}/run_inference.py\
            inference.output_prefix={rf_dir}/pred\
            inference.input_pdb={input_pdb}\
            'contigmap.contigs=[A1-{rid_init-1}/{rid_fin - rid_init}-{rid_fin - rid_init}/A{rid_fin}-{len_heavy}/0 B1-{len_light}]'\
            diffuser.partial_T={partial_steps}\
            'contigmap.provide_seq=[0-{rid_init-1},{rid_fin}-{len_heavy},{len_heavy+1}-{len_heavy+len_light}]'\
            inference.num_designs={num_designs}\
            """
    os.system(command)

def rfdiff_full(input_pdb,
                rf_dir,
                rid_init,
                rid_fin,
                len_heavy,
                len_light,
                len_ant,
                rf_script_path,
                se3_env,
                partial_steps=10,
                num_designs=1,
                ):

    command = f"""module use /soft/modulefiles &&\
            module load conda &&\
            conda activate {se3_env}\
            && {rf_script_path}/run_inference.py\
            inference.output_prefix={rf_dir}/pred\
            inference.input_pdb={input_pdb}\
            'contigmap.contigs=[A1-{rid_init-1}/{rid_fin - rid_init}-{rid_fin - rid_init}/A{rid_fin}-{len_heavy}/0 B1-{len_light}/0 C1-{len_ant}]'\
            diffuser.partial_T={partial_steps}\
            'contigmap.provide_seq=[0-{rid_init-1},{rid_fin}-{len_heavy},{len_heavy+1}-{len_heavy+len_light},{len_heavy+len_light+1}-{len_heavy+len_light+len_ant}']\
            inference.num_designs={num_designs}\
            """
    os.system(command)
