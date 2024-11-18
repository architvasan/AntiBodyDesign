#!/bin/bash
#PBS -N rfdiff_14
#PBS -l select=1
#PBS -q by-node
#PBS -l walltime=10:00:00
#PBS -A datascience
#PBS -l filesystems=eagle
#PBS -m abe
#PBS -M avasan@anl.gov

module use /soft/modulefiles
module load conda

conda activate /lus/eagle/projects/datascience/avasan/envs/chai1
conda env list
#module load mpiwrappers/cray-mpich-llvm 
#module load PrgEnv-gnu 
#export CRAY_ACCEL_TARGET="nvidia80" 
#export CRAY_TCMALLOC_MEMFS_FORCE="1" 
#export CRAYPE_LINK_TYPE="dynamic" 
#export CRAY_ACCEL_VENDOR="nvidia"
#export CRAY_CPU_TARGET="x86-64"

trial=2
cd /lus/eagle/projects/datascience/avasan/Simulations/AntiBodyDesign

stdout=logs/run_t${trial}.mpi.a3HFM_b4NCO.afterdpo_norf.log
stderr=logs/run_t${trial}.mpi.a3HFM_b4NCO.afterdpo_norf.err
NDEPTH=16

data_dir=iedb_tables_pos_and_neg
output_dir_general=trials/T${trial}_ant3HFM_body4NCO_afterdpo_norf
log_name=logs/run_t${trial}_a3HFM_b4NCO_afterdpo_norf
dbfil=fab_lyso.csv
chaintarget=heavy_chain
cdrid=heavy_cdr3
ndesigns=20

totalrank=8
pernode=8

mkdir $output_dir_general

# Run each rank separately on sophia since there are issues with
# mpi4py

for i in $(seq 0 6);
do
    rank=${i}
    echo $i
    /lus/eagle/projects/datascience/avasan/envs/chai1/bin/python run_pipeline_full_newcdrs_nordfiff.py \
        --dbfil $data_dir/$dbfil \
        --cdrfil $data_dir/generated_seqs_dpo_${rank}_rst0.csv \
        --chaintarget $chaintarget \
        --cdrid $cdrid \
        --outdir $output_dir_general/${rank} \
        --ndesigns $ndesigns \
        -G $pernode \
        -D ${rank} > $log_name.${rank}.log 2> $log_name.$rank.err &
done    

rank=7
echo $i
    /lus/eagle/projects/datascience/avasan/envs/chai1/bin/python run_pipeline_full_newcdrs_nordfiff.py \
        --dbfil $data_dir/$dbfil \
        --cdrfil $data_dir/generated_seqs_dpo_${rank}_rst0.csv \
        --chaintarget $chaintarget \
        --cdrid $cdrid \
        --outdir $output_dir_general/${rank} \
        --ndesigns $ndesigns \
        -G $pernode \
        -D ${rank} > $log_name.${rank}.log 2> $log_name.$rank.err 




#mpiexec -n $totalrank -ppn $pernode \
#    --depth=${NDEPTH} --cpu-bind depth \
#    ./set_affinity_gpu_polaris.sh \
#    /lus/eagle/projects/datascience/avasan/envs/chai1/bin/python run_pipeline_full_newcdrs.py \
#    --dbfil $data_dir/$dbfil \
#    --cdrfil $data_dir/heavy-cdr3_all.json \
#    --chaintarget $chaintarget \
#    --usempi \
#    --cdrid $cdrid \
#    --outdir $output_dir_general \
#    --ndesigns $ndesigns \
#    -G $pernode \
#    > $stdout 2> $stderr
