#!/bin/bash

#SBATCH --job-name=cellranger_count
#SBATCH --array=1-6
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 72
#SBATCH -t 2:00:00
#SBATCH -o ./logs/cellranger_"stdout.%j.%N"
#SBATCH -e ./logs/cellranger_"stderr.%j.%N"
#SBATCH --mem=100G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=

set -e
mkdir -p cellranger_out
# mkdir -p cellranger_out_10000

# this should create an array of jobs that performs all 6 individual mappings

echo $SLURM_ARRAY_TASK_ID

# this returns the line from the bio_samples file that corresponds to the 
# integer from the $SLURM_ARRAY_TASK_ID variable (which is an integer 1-6)
SAMPLE_ID=$(cat outputs/bio_samples.txt | sed -n ${SLURM_ARRAY_TASK_ID}p) 

echo $SAMPLE_ID


source ~/.bashrc

cellranger count --id="$SAMPLE_ID" \
                --transcriptome=Bos_taurus.ARS-UCD1.2 \
                --fastqs=raw_data \
                --sample="$SAMPLE_ID"-1,"$SAMPLE_ID"-2,"$SAMPLE_ID"-3,"$SAMPLE_ID"-4 \
                --localcores=72 \
                --disable-ui

mv "$SAMPLE_ID"/ ./cellranger_out/

#####

# This block is what I would use to try to capture neutrophils as reccommended by
# cell ranger documentation.  However, downstream scripts will need to be modified
# to distinguish empty cells from those with a low number of transcripts detected.

# cellranger count --id="$SAMPLE_ID" \
#                  --transcriptome=Bos_taurus.ARS-UCD1.2 \
#                  --fastqs=raw_data \
#                  --sample="$SAMPLE_ID"-1,"$SAMPLE_ID"-2,"$SAMPLE_ID"-3,"$SAMPLE_ID"-4 \
#                  --localcores=72 \
#                  --disable-ui \
# 		 --force-cells=10000
#
# mv "$SAMPLE_ID"/ ./cellranger_out_10000/

#####


# start the R scripts
#sbatch scripts/run_R_scripts.slurm

