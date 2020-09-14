#!/bin/bash

# Submit this script with: sbatch <this-filename>
#SBATCH --account=seinfeldgroup 		# accounting resource
#SBATCH --job-name=test_varN			# job name

#SBATCH --time=00:30:00				# walltime
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
MU1=0.05e-6
SIG1=1.2
N1=50e6
KAP1=0.7

MU2=0.01e-6
SIG2=1.2
N2=0.0
KAP2=0.7

# run time params
SPIN=300
NT=600
OUTFREQ=300

# case params
NX=129
NZ=301

# output directory
echo array job id = ${SLURM_ARRAY_JOB_ID}
echo array task id = ${SLURM_ARRAY_TASK_ID}
OUTPUT_DIR=$HOME/output/test_arrays/mu_${MU1}_kappa_${KAP1}_N_${N1}/exp${SLURM_ARRAY_TASK_ID}
#OUTPUT_DIR=/central/scratch/csinger/tmp_output_lgr/test_arrays/mu_${MU1}_kappa_${KAP1}_N_${N1}/exp${SLURM_ARRAY_TASK_ID}

# set param strings
CASE_PARAMS="--case=dycoms_rf02 --nx=$NX --ny=0 --nz=$NZ"
RUN_PARAMS="--dt=1 --spinup=$SPIN --nt=$NT --outfreq=$OUTFREQ --backend=CUDA"
COMMON_PARAMS="--sstp_cond=10 --sstp_coal=10 --rng_seed=0"
MICRO_PARAMS="--micro=lgrngn --mean_rd1=$MU1 --sdev_rd1=$SIG1 --n1_stp=$N1 --kappa1=$KAP1 --mean_rd2=$MU2 --sdev_rd2=$SIG2 --n2_stp=$N2 --kappa2=$KAP2 --sd_conc=40"

wet_bins_str=$(python make_bins.py "wet")
dry_bins_str=$(python make_bins.py "dry")

BINS="--out_wet=$wet_bins_str --out_dry=$dry_bins_str"

# copy this input script to the output directory
mkdir -p $OUTPUT_DIR
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
mv slurm-${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out $OUTPUT_DIR/slurm.out
