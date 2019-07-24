#!/bin/bash

# Submit this script with: sbatch <this-filename>
#SBATCH --account=seinfeldgroup 		# accounting resource
#SBATCH --job-name=test-dycoms			# job name
#SBATCH --output=slurm.out			# output filename

#SBATCH --time=00:10:00				# walltime
#SBATCH --nodes=1				# 1 GPU node
#SBATCH --ntasks=1				# 1 CPU core to drive GPU
#SBATCH --gres=gpu:1				# request 1 GPU

##SBATCH --exclusive				# for exclusive use of node
##SBATCH --qos=debug				# debug mode

#SBATCH --mail-user=csinger@caltech.edu   	# email address
#SBATCH --mail-type=ALL				# notify all

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
module load singularity/3.2.1
module list

# pathes to directories
BASE=/home/csinger/microphys
CONTAINER=$BASE/sng_output.sif
MODEL=$BASE/UWLCM/build/src/bicycles

# aerosol distribution params
MU1=0.011e-6
SIG1=1.2
N1=125e6
KAP1=0.61

MU2=0.06e-6
SIG2=1.7
N2=65e6
#KAP2=0.61

# run time params
SPIN=1800
NT=3600
OUTFREQ=300

# case params
NX=33 #129
NZ=76 #301

# output directory
OUTPUT_DIR=$BASE/output_lgr/dycoms/spin_${SPIN}_nt_${NT}_out_${OUTFREQ}

# set param strings
CASE_PARAMS="--case=dycoms_rf02 --nx=$NX --ny=0 --nz=$NZ"
RUN_PARAMS="--dt=1 --spinup=$SPIN --nt=$NT --outfreq=$OUTFREQ --backend=CUDA"
COMMON_PARAMS="--sstp_cond=10 --sstp_coal=10 --rng_seed=0"
MICRO_PARAMS="--micro=lgrngn --mean_rd1=$MU1 --sdev_rd1=$SIG1 --n1_stp=$N1 --kappa1=$KAP1 --mean_rd2=$MU2 --sdev_rd2=$SIG2 --n2_stp=$N2 --kappa2=$KAP1 --sd_conc=40"

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
srun singularity -d exec --nv $CONTAINER $MODEL --outdir=$OUTPUT_DIR $CASE_PARAMS $RUN_PARAMS $COMMON_PARAMS $MICRO_PARAMS $BINS

# copy slurm file to output directory
mv slurm.out $OUTPUT_DIR/slurm.out
