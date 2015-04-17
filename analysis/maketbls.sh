#!/bin/sh

timestamp=$1
start=$2
inc=$3
end=$4

for index in `seq $start $inc $end`; do

echo index = $index

results_path=`printf "DGC_forest_%d_test_results_passive_$timestamp.dat" $index`
results_dir=/som/iraikov/DGC_forest_test_results/passive/$timestamp/$index

echo results_path=$results_path
echo results_dir=$results_dir

echo ../../analysis/maketbl.sh $results_path $results_dir
../../analysis/maketbl.sh $results_path $results_dir

done

