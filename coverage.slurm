#!/bin/bash

# check that the script is launched with sbatch
if [ "x$SLURM_JOB_ID" == "x" ]; then
   echo "You need to submit your job to the queuing system with sbatch"
   exit 1
fi

# Partition for the job:
#SBATCH -p main

# The name of the job:
#SBATCH --job-name="coverage"

# Maximum number of CPU cores used by the job:
#SBATCH --ntasks=1

# The amount of memory in megabytes per process in the job:
#SBATCH --mem-per-cpu=8192

# Send yourself an email when the job:
# aborts abnormally (fails)
#SBATCH --mail-type=FAIL
# begins
#SBATCH --mail-type=BEGIN
# ends successfully
#SBATCH --mail-type=END

# The maximum running time of the job in days-hours:mins:sec
#SBATCH --time=0-10:00:00

# Run the job from the directory where it was launched (default):
# The modules to load:
module load bedtools-intel/2.22.1 

declare -a BAMS=( $XXX_PUT_LIST_OF_BAMS_IN_HERE ); 

BAMIN=${BAMS[$SLURM_ARRAY_TASK_ID]}
OUT=${BAMIN}.coverage

bedtools coverage -d -abam $BAMIN -b final_gene_list.bed > $OUT
