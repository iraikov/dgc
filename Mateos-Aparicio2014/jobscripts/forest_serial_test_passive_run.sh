#!/bin/bash
#
#$ -q som,asom,free64,pub64
#$ -t 1-128
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -N DGC_forest_test_MateosAparicio2014_passive_na8st
#$ -o ./results/forest_test_MateosAparicio2014_passive_na8st.$JOB_ID.o


module load neuron/7.4alpha

set -x

results_path="/pub/iraikov/DGC_forest_test_results/`date +'%Y%m%d%H%M'`"
export results_path
mkdir -p $results_path

nrniv -nobanner -nogui -c "batch_size=100" -c "task_id=$SGE_TASK_ID - 1" \
-c "strdef forest_config" -c 'forest_config="./config/DGC_forest_hpc_passive.config"' \
./DGC_serial_test_from_forest_passive_na8st.hoc
