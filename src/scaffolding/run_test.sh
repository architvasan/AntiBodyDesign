#!/bin/bash
# Check if the correct number of arguments is provided

module use /soft/modulefiles
module load conda
conda activate /lus/eagle/projects/datascience/avasan/envs/SE3nv

src_dir=/eagle/datascience/avasan/RFDiffusionProject/RFdiffusion/scripts/
# Try a command
#if mkdir logs/; then
#  echo "made logs"
#else
#  echo "logs exists"
#fi

#CUDA_VISIBLE_DEVICES=$device $src_dir/run_inference.py inference.output_prefix=${output_prefix} inference.input_pdb=${input_pdb} "contigmap.contigs=[A1-${size_ref}/0 ${min_design}-${max_design}]" "ppi.hotspot_res=[${hotspot_res}]" inference.num_designs=$num_designs denoiser.noise_scale_ca=0 denoiser.noise_scale_frame=0 > logs/${outfile}.log 2> logs/${outfile}.err

CUDA_VISIBLE_DEVICES=0 $src_dir/run_inference.py inference.output_prefix=./tests/test_rf.pdb inference.input_pdb=./tests/pred.model_idx_0.pdb 'contigmap.contigs=[10-40/A1-96/A116-227]' inference.num_designs=10
#../scripts/run_inference.py inference.output_prefix=example_outputs/design_partialdiffusion_peptidewithsequence inference.input_pdb=input_pdbs/peptide_complex_ideal_helix.pdb 'contigmap.contigs=["172-172/0 34-34"]' diffuser.partial_T=10 inference.num_designs=10 'contigmap.provide_seq=[172-205]'
