#!/bin/bash

module use /soft/modulefiles
module load conda

cd /eagle/datascience/avasan/Simulations/AntiBodyDesign
conda activate /lus/eagle/projects/datascience/avasan/envs/chai1
conda env list
#module load mpiwrappers/cray-mpich-llvm 

inputfil=iedb_tables_pos_and_neg/hiv_ab.csv
chaiout=trials/T1/chaiout
rfout=trials/T2/rfout
logout=trials/T2/logout
resmaploc=trials/T1/resmaps
numdesign=10
foldinit=False
totalrank=4
pernode=4
stdout=logs/rf_t2.mpi.sophia.log
stderr=logs/rf_t2.mpi.sophia.err


mpiexec -n $totalrank -ppn $pernode \
    --depth=8 \
    --cpu-bind depth \
    ./set_affinity_gpu_polaris.sh \
    python run_pipeline_full_mpi4py.py \
        -i $inputfil \
        -C $chaiout \
        -R $rfout \
        -L $logout \
        -M $resmaploc \
        -N $numdesign \
        -G $pernode \
        -F $foldinit \
            > $stdout 2> $stderr