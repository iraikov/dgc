#!/bin/bash
#SBATCH -J DGC_forest_synapse_test_MateosAparicio2014_na8st
#SBATCH -A TG-IBN140007
#SBATCH -o ./results/DGC_forest_synapse_test_MateosAparicio2014_na8st.%j.o
#SBATCH -N 1 -n 1
#SBATCH -p serial
#SBATCH -t 01:00:00
#SBATCH --mail-user=ivan.g.raikov@gmail.com
#SBATCH --mail-type=END

set -x
echo TASK_ID=$TASK_ID
x86_64/special -mpi -c "batch_size=6" -c "task_id=$TASK_ID" \
-c "strdef forest_config" -c 'forest_config="./config/DGC_forest_hpc.config"' \
./DGC_serial_synapse_test_from_forest_na8st.hoc
