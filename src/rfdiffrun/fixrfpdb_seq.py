import os

def fixpdb(rf_dir,
           dlbind_env,
           helpersrc,
           ):
    
    command = f"""\
              module use /soft/modulefiles &&\
              module load conda &&\
              conda activate {dlbind_env} &&\
              python {helpersrc}/addFIXEDlabels.py\
              --pdbdir  {rf_dir}\
              --trbdir {rf_dir}\
              --verbose\
              """
    os.system(command)
