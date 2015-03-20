#!/bin/bash
#
#$ -q som,asom,free64,pub64
#$ -t 1-20
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -N DGC_forest_synapse_test_MateosAparicio2014_na8st
#$ -o ./results/forest_synapse_test_MateosAparicio2014_na8st.$JOB_ID.o

module load neuron/7.4alpha
results_path="/pub/iraikov/DGC_forest_synapse_test_results/`date +'%Y%m%d%H%M'`"
export results_path
echo nrniv -nobanner -nogui -c "batch_size=20" -c "task_id=$SGE_TASK_ID - 1" \
-c "strdef forest_config" -c 'forest_config="./config/DGC_forest_hpc.config"' \
./DGC_serial_synapse_test_from_forest_na8st.hoc

nrniv -nobanner -nogui -c "batch_size=20" -c "task_id=$SGE_TASK_ID - 1" \
-c "strdef forest_config" -c 'forest_config="./config/DGC_forest_hpc.config"' \
./DGC_serial_synapse_test_from_forest_na8st.hoc
