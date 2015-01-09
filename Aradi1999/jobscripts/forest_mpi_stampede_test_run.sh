#!/bin/bash

#SBATCH -J DGC_forest_test # name of job
#SBATCH -o DGC_forest_test.%j.o # output file (errors can be placed in this file, too). %j expands to the job id number.
#SBATCH -n 20 # number of cores
#SBATCH -p normal # queue to use
#SBATCH -t 02:00:00 # time limit
#SBATCH --mail-user=ivan.g.raikov@gmail.com # where to send notifications
#SBATCH --mail-type=ALL # types of mail notifications to receive

ibrun tacc_affinity ./x86_64/special -mpi - <<EOF
strdef forest_config
forest_config="DGC_stampede_forest_config.hoc"
load_file("DGC_mpi_test_from_forest_na8st.hoc")
EOF

