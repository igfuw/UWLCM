#!/bin/bash

# Submit this script with: sbatch <this-filename>
#SBATCH --account=seinfeldgroup 		# accounting resource
#SBATCH --job-name=run-dycoms			# job name

#SBATCH --time=05:00:00				# walltime
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
MU1=0.03e-6
SIG1=1.28
N1=90e6
KAP1=0.61

MU2=0.14e-6
SIG2=1.75
N2=15e6
KAP2=0.61

# run time params
SPIN=3600
NT=86400
OUTFREQ=600

# case params
NX=129
NZ=101

# output directory
OUTPUT_DIR=$HOME/output/rico/exp1

# set param strings
CASE_PARAMS="--case=rico --nx=$NX --ny=0 --nz=$NZ"
RUN_PARAMS="--dt=1 --spinup=$SPIN --nt=$NT --outfreq=$OUTFREQ --backend=CUDA"
COMMON_PARAMS="--sstp_cond=10 --sstp_coal=10 --rng_seed=0"
MICRO_PARAMS="--micro=lgrngn --mean_rd1=$MU1 --sdev_rd1=$SIG1 --n1_stp=$N1 --kappa1=$KAP1 --mean_rd2=$MU2 --sdev_rd2=$SIG2 --n2_stp=$N2 --kappa2=$KAP2 --sd_conc=40"

wet_bins_str=$(python make_bins.py "wet")
dry_bins_str=$(python make_bins.py "dry")

BINS="--out_wet=$wet_bins_str --out_dry=$dry_bins_str"

# copy this input script to the output directory
mkdir $OUTPUT_DIR
cp ${0} $OUTPUT_DIR/input.sh

# make .txt file in output directory with wet_bins_str and dry_bins_str
BINS_FILE="bin_strs.txt"
echo wet_bins_str: >> $OUTPUT_DIR/$BINS_FILE
echo $wet_bins_str >> $OUTPUT_DIR/$BINS_FILE
echo dry_bins_str: >> $OUTPUT_DIR/$BINS_FILE
echo $dry_bins_str >> $OUTPUT_DIR/$BINS_FILE
chmod -w $OUTPUT_DIR/$BINS_FILE

# run program (-d == in debug mode)
echo about to call bicycles...
echo job id = ${SLURM_JOB_ID}
srun singularity exec --nv $CONTAINER $MODEL --outdir=$OUTPUT_DIR $CASE_PARAMS $RUN_PARAMS $COMMON_PARAMS $MICRO_PARAMS $BINS
echo done!

# copy slurm file to output directory
mv slurm-${SLURM_JOB_ID}.out $OUTPUT_DIR/slurm.out
