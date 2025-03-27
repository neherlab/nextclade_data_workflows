#!/usr/bin/env bash

#SBATCH --output=log/%j.out                 # where to store the output ( %j is the JOBID )
#SBATCH --error=log/%j.err                  # where to store error messages

# activate conda environment
conda activate py310simple
. /scicore/home/neher/roemer0001/miniconda3/etc/profile.d/conda.sh
export AUGUR_MINIFY_JSON=1
export AUGUR_RECURSION_LIMIT=10000

{exec_job}
