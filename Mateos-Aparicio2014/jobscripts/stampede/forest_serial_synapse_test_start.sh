#!/bin/bash

config_path=./config/DGC_forest_hpc.config

ncells=1000
data_path=$WORK/DGC_forest/dat
results_path="./results/DGC_forest_synapse_test_results/`date +'%Y%m%d%H%M'`"

echo "$ncells" > $config_path
echo "$data_path" >> $config_path
echo "$results_path" >> $config_path

for i in `seq 1 6`; do
  TASK_ID=$i
  export TASK_ID
  sbatch ./jobscripts/stampede/forest_serial_synapse_test_run.sh
done


