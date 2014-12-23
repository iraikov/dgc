#!/bin/bash
#
#$ -q som,asom,pub64,free64
#$ -pe openmp 64
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -N DGC_forest_test
#$ -o ./results/forest_test_run.$JOB_ID.o
#$ -R y

module load neuron/7.3
mpiexec -np 64 nrniv -mpi -nobanner -nogui ./DGC_test_from_forest_na8st.hoc


