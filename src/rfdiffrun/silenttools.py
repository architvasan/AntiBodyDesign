import os

def pdb2silent(dlbind_env,
               rf_dir,
               silent_tool_path,
                ):

    command = f"""\
             module use /soft/modulefiles &&\
             module load conda &&\
             source ~/.bashrc &&\
             conda activate {dlbind_env} &&\
             {silent_tool_path}/silentfrompdbs\
             {rf_dir}/*.pdb >\
             {rf_dir}/rfout.silent\
             """
    
    os.system(command)

def extractpdb(localpath,
               dlbind_env,
               rf_dir,
               silent_tool_path):
    command = f"""\
              module use /soft/modulefiles &&\
              module load conda &&\
              source ~/.bashrc &&\
              cd {rf_dir} &&\
              conda activate {dlbind_env} &&\
              {silent_tool_path}/silentextract\
              mpnnout.silent &&\
              cd /eagle/datascience/avasan/Simulations/Antibody_Design\
              """
    os.system(command)


    return 
