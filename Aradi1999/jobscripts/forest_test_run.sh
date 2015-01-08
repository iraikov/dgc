#!/bin/bash
#
#$ -q free64,pub64
#$ -t 1-20
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -N DGC_forest_test_Aradi1999_na8st
#$ -o ./results/forest_test_Aradi1999_na8st.$JOB_ID.o


module load neuron/7.3
nrniv -nobanner -nogui -c "batch_size=20" -c "task_id=$SGE_TASK_ID - 1" ./DGC_test_from_forest_na8st.hoc


