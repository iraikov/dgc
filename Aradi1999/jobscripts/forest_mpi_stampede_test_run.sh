#!/bin/bash

#SBATCH -J DGC_forest_test # name of job
#SBATCH -o DGC_forest_test.%j.o # output file (errors can be placed in this file, too). %j expands to the job id number.
#SBATCH -n 1 # number of cores
#SBATCH -p development # queue to use
#SBATCH -t 00:10:00 # time limit
#SBATCH --mail-user=ivan.g.raikov@gmail.com # where to send notifications
#SBATCH --mail-type=ALL # types of mail notifications to receive

ibrun tacc_affinity ./x86_64/special -mpi -c "strdef forest_config" -c "forest_config=\"DGC_stampede_forest_config.hoc\"" \
 ./DGC_mpi_test_from_forest_na8st.hoc
