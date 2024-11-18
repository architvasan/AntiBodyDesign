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
#module load PrgEnv-gnu 
#export CRAY_ACCEL_TARGET="nvidia80" 
#export CRAY_TCMALLOC_MEMFS_FORCE="1" 
#export CRAYPE_LINK_TYPE="dynamic" 
#export CRAY_ACCEL_VENDOR="nvidia"
#export CRAY_CPU_TARGET="x86-64"

trial=16
cd /eagle/datascience/avasan/Simulations/AntiBodyDesign
NDEPTH=16

log_name=logs/run_t${trial}_a3HFM_b4NCO.sim

totalrank=8
pernode=8

cd /eagle/datascience/avasan/Simulations/AntiBodyDesign_Constrained
# Run each rank separately on sophia since there are issues with
# mpi4py

for i in $(seq 0 6);
do
    rank=${i}
    echo $i
    mkdir trials/T1_ant3HFM_body4NCO_simulations/$i
    /lus/eagle/projects/datascience/avasan/envs/chai1/bin/python src/simulate/simulate_implicit.py \
        -df allpdbs_t16_${i}.csv \
        -o trials/T16_ant3HFM_body4NCO_simulations/$i \
        -t 250000 \
        -D ${rank} > $log_name.${rank}.log 2> $log_name.$rank.err &
done    

rank=7
echo $i
i=7

mkdir trials/T16_ant3HFM_body4NCO_simulations/$i
/lus/eagle/projects/datascience/avasan/envs/chai1/bin/python src/simulate/simulate_implicit.py \
    -df allpdbs_t16_${i}.csv \
    -o trials/T16_ant3HFM_body4NCO_simulations/$i \
    -t 250000 \
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
