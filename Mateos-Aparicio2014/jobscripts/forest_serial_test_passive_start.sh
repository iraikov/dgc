#!/bin/bash

config_path=./config/DGC_forest_hpc_passive.config

ncells=1000
data_path_prefix=/som/iraikov/DGC_forest
results_path_prefix="/som/iraikov/DGC_forest_test_results/passive/`date +'%Y%m%d%H%M'`"
forest_seq="10 10 10"

echo "$ncells" > $config_path
echo "$data_path_prefix" >> $config_path
echo "$results_path_prefix" >> $config_path
echo "$forest_seq" >> $config_path

qsub ./jobscripts/forest_serial_test_passive_run.sh
