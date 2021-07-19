#!/bin/bash

# set job id
JOBNO=032

# clean run directory and link all required files
./prep_run_a2.sh

# submit the job
sbatch -J AISOMIP_$JOBNO \
       --account $HECACC \
       run_a2.sh

