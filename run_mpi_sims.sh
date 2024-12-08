#!/bin/bash
#PBS -N rfdiff_14
#PBS -l select=10
#PBS -l walltime=1:00:00
#PBS -q debug-scaling
#PBS -A datascience
#PBS -l filesystems=eagle
#PBS -m abe
#PBS -M avasan@anl.gov

module use /soft/modulefiles
module load conda

conda activate /lus/eagle/projects/datascience/avasan/envs/chai1
conda env list
#module load mpiwrappers/cray-mpich-llvm 
module load PrgEnv-gnu 
export CRAY_ACCEL_TARGET="nvidia80" 
export CRAY_TCMALLOC_MEMFS_FORCE="1" 
export CRAYPE_LINK_TYPE="dynamic" 
export CRAY_ACCEL_VENDOR="nvidia"
export CRAY_CPU_TARGET="x86-64"

trial=14
cd /eagle/datascience/avasan/Simulations/AntiBodyDesign
stdout=logs/sims_9_10_14.log
stderr=logs/sims_9_10_14.err
NDEPTH=16

data_dir=iedb_tables_pos_and_neg
output_dir_general=trials/T${trial}_ant3HFM_body4NCO
log_name=logs/run_t${trial}_a3HFM_b4NCO
dbfil=fab_lyso.csv
chaintarget=heavy_chain
cdrid=heavy_cdr3
ndesigns=10

totalrank=40
pernode=4

#mkdir $output_dir_general

mpirun -n $totalrank -ppn $pernode \
    --depth=${NDEPTH} --cpu-bind depth \
    /lus/eagle/projects/datascience/avasan/envs/chai1/bin/python simulate_structures.py \
    -df all_pdbs_save_loc.csv -o trials/T9_10_14_antdy4NCO_simulations/ --usempi \
    > $stdout 2> $stderr


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
#    > $stdout 2> $stderr
