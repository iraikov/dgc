#!/bin/bash

ncells=1000
data_path=/data/users/iraikov/model/DGC_forest/dat
results_path="/pub/iraikov/DGC_forest_test_results/`date +'%Y%m%d%H%M'`"

echo "$ncells" > DGC_forest_hpc.config
echo "$data_path" >>  DGC_forest_hpc.config
echo "$results_path" >>  DGC_forest_hpc.config

qsub ./jobscripts/forest_test_run.sh
