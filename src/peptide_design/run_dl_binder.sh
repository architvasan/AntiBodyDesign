#!/bin/bash
# Check if the correct number of arguments is provided
if [ "$#" -lt 8 ]; then
  echo "Usage: $0 path_mpnn path_af2 input_silent mpnn_out af2_out device checkpoint log"
  exit 1
fi

path_mpnn=$1
path_af2=$2
input_silent=$3
mpnn_out=$4
af2_out=$5
device=$6
checkpoint=$7
log=$8

module use /soft/modulefiles
module load conda
conda activate /lus/eagle/projects/datascience/avasan/envs/proteinmpnn_binder_design

#echo "${size_ref}"

CUDA_VISIBLE_DEVICES=$device ${path_mpnn}/dl_interface_design.py -silent $input_silent -outsilent $mpnn_out -checkpoint_name ${checkpoint}_mpnn > $log.log 2> $log.err

##################################### AlphaFold2 Complex Prediction ##########################################
##############################################################################################################

conda deactivate
conda activate /lus/eagle/projects/datascience/avasan/envs/af2_binder_design
CUDA_VISIBLE_DEVICES=$device ${path_af2}/predict.py -silent $mpnn_out -outsilent $af2_out -checkpoint_name $checkpoint > $log.af2.log 2> $log.af2.err
