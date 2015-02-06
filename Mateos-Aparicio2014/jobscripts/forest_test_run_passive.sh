#!/bin/bash
#
#$ -q free64,pub64
#$ -t 1-20
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -N DGC_forest_test_MateosAparicio2014_passive_na8st
#$ -o ./results/forest_test_MateosAparicio2014_passive_na8st.$JOB_ID.o


module load neuron/7.3
mkdir -p ./results
nrniv -nobanner -nogui -c "batch_size=20" -c "task_id=$SGE_TASK_ID - 1" \
-c "strdef forest_config" -c 'forest_config="DGC_hpc_forest_config_passive.hoc"' \
./DGC_serial_test_from_forest_passive_na8st.hoc
