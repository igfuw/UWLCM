#!/bin/bash

# Submit this script with: sbatch <this-filename>
#SBATCH --account=seinfeldgroup 		# accounting resource
#SBATCH --job-name=test-blk2			# job name

#SBATCH --time=01:00:00				# walltime
#SBATCH --nodes=1				# 1 GPU node
#SBATCH --ntasks=1				# 1 CPU core to drive GPU
#SBATCH --gres=gpu:1				# request 1 GPU

#SBATCH --mail-user=csinger@caltech.edu   	# email address
#SBATCH --mail-type=ALL				# notify all

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
module load singularity/3.2.1
module list

# pathes to directories
HOME=/home/csinger/microphys
CONTAINER=$HOME/sng_output.sif
MODEL=$HOME/UWLCM/build/src/bicycles

# aerosol distribution params
# dycoms
MU1=0.01e-6
SIG1=1.2
N1=100e6
KAP1=1.3

# run time params
SPIN=3600
NT=21600
OUTFREQ=600

# case params
NX=129
NZ=301

# output directory
OUTPUT_DIR=$HOME/output/test_blk_2m/mu_${MU1}_kappa_${KAP1}_N_${N1}/exp42

# set param strings
CASE_PARAMS="--case=dycoms_rf02 --nx=$NX --ny=0 --nz=$NZ"
RUN_PARAMS="--dt=1 --spinup=$SPIN --nt=$NT --outfreq=$OUTFREQ --backend=CUDA"
COMMON_PARAMS="--sstp_cond=10 --sstp_coal=10 --rng_seed=42"
TURB_PARAMS="--sgs=0 --turb_adve=0"
MICRO_PARAMS="--micro=blk_2m --mean_rd1=$MU1 --sdev_rd1=$SIG1 --n1_stp=$N1 --kappa1=$KAP1"

#wet_bins_str=$(python make_bins.py "wet")
#dry_bins_str=$(python make_bins.py "dry")
#
#BINS="--out_wet=$wet_bins_str --out_dry=$dry_bins_str"
#echo $BINS

# copy this input script to the output directory
mkdir -p $OUTPUT_DIR
cp ${0} $OUTPUT_DIR/input.sh

## make .txt file in output directory with wet_bins_str and dry_bins_str
#BINS_FILE="bin_strs.txt"
#echo wet_bins_str: >> $OUTPUT_DIR/$BINS_FILE
#echo $wet_bins_str >> $OUTPUT_DIR/$BINS_FILE
#echo dry_bins_str: >> $OUTPUT_DIR/$BINS_FILE
#echo $dry_bins_str >> $OUTPUT_DIR/$BINS_FILE
#chmod -w $OUTPUT_DIR/$BINS_FILE

# run program (-d == in debug mode)
echo about to call bicycles...
echo job id = ${SLURM_JOB_ID}
srun singularity exec --nv $CONTAINER $MODEL --outdir=$OUTPUT_DIR $CASE_PARAMS $RUN_PARAMS $COMMON_PARAMS $TURB_PARAMS $MICRO_PARAMS $BINS
echo done!

# copy slurm file to output directory
mv slurm-${SLURM_JOB_ID}.out $OUTPUT_DIR/slurm.out
