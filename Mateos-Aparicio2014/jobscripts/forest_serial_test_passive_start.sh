#!/bin/bash

config_path=./config/DGC_forest_hpc_passive.config

ncells=1000
data_path=/data/users/iraikov/model/DGC_forest/dat
results_path="/pub/iraikov/DGC_forest_test_results/passive/`date +'%Y%m%d%H%M'`"

echo "$ncells" > $config_path
echo "$data_path" >> $config_path
echo "$results_path" >> $config_path

qsub ./jobscripts/forest_serial_test_passive_run.sh
