#!/bin/sh

#SBATCH --output=log/%j.out                 # where to store the output ( %j is the JOBID )
#SBATCH --error=log/%j.err                  # where to store error messages

# activate conda environment
source /scicore/home/neher/roemer0001/miniconda3/etc/profile.d/conda.sh
conda activate nextstrain

{exec_job}
