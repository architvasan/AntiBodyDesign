#!/bin/bash
# Check if the correct number of arguments is provided
if [ "$#" -lt 6 ]; then
  echo "Usage: $0 path_af2 mpnn_out af2_out device checkpoint log"
  exit 1
fi

path_af2=$1
mpnn_out=$2
af2_out=$3
device=$4
checkpoint=$5
log=$6

module use /soft/modulefiles
module load conda
##################################### AlphaFold2 Complex Prediction ##########################################
##############################################################################################################

conda deactivate
conda activate /lus/eagle/projects/datascience/avasan/envs/af2_binder_design
echo ${path_af2}
CUDA_VISIBLE_DEVICES=$device ${path_af2}/predict.py -silent $mpnn_out -outsilent $af2_out -checkpoint_name $checkpoint > $log.af2.log 2> $log.af2.err 
