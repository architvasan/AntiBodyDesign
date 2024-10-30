#!/bin/bash
pdb_dir=$1
out_silent=$2
silentpwd=/lus/eagle/projects/datascience/avasan/RFDiffusionProject/dl_binder_design/include/silent_tools

module use /soft/modulefiles
module load conda
conda activate /lus/eagle/projects/datascience/avasan/envs/proteinmpnn_binder_design
# pdbs to silent
${silentpwd}/silentfrompdbs ${pdb_dir}/*.pdb > ${out_silent}
