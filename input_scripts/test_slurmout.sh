#!/bin/bash

# Submit this script with: sbatch <this-filename>
#SBATCH --account=seinfeldgroup                 # accounting resource
#SBATCH --job-name=test	                 	# job name
##SBATCH --output=slurm.out			# output filename

#SBATCH --time=00:01:00                         # walltime
#SBATCH --nodes=1                               # 1 GPU node
#SBATCH --ntasks=1                              # 1 CPU core to drive GPU
#SBATCH --gres=gpu:1                            # request 1 GPU

##SBATCH --exclusive                            # for exclusive use of node
##SBATCH --qos=debug                            # debug mode

#SBATCH --mail-user=csinger@caltech.edu         # email address
#SBATCH --mail-type=ALL                         # notify all

BASE=/home/csinger/microphys/
OUTPUT_DIR=$BASE/output_lgr/dycoms/test_slurm_copy
mkdir $OUTPUT_DIR

echo hello
echo goodbye
echo $SLURM_JOB_ID

mv slurm-${SLURM_JOB_ID}.out $OUTPUT_DIR/slurm.out
