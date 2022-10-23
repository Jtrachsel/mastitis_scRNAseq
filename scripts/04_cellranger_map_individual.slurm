#!/bin/bash

#SBATCH --job-name=cellranger_count
#SBATCH --array=1-6                              # name of the job submitted
#SBATCH -p short                                        # name of the queue you are submitting to
#SBATCH -N 1                                            # number of nodes in this job
#SBATCH -n 72                                           # number of cores/tasks in this job, you get all 20 cores with 2 threads per core with hyperthreading
#SBATCH -t 2:00:00                                     # time allocated for this job hours:mins:seconds
#SBATCH -o ./logs/cellranger_"stdout.%j.%N"                               # standard out %j adds job number to outputfile name and %N adds the node name
#SBATCH -e ./logs/cellranger_"stderr.%j.%N"                               # optional but it prints our standard error
#SBATCH --mem=100G   
#SBATCH --mail-type=ALL
#SBATCH --mail-user=julestrachsel@gmail.com

set -e
mkdir -p cellranger_out
mkdir -p cellranger_out_10000

# this should create an array of jobs that performs all 6 individual mappings

echo $SLURM_ARRAY_TASK_ID

SAMPLE_ID=$(cat outputs/bio_samples.txt | sed -n ${SLURM_ARRAY_TASK_ID}p)

echo $SAMPLE_ID


source ~/.bashrc

#cellranger count --id="$SAMPLE_ID" \
#                 --transcriptome=Bos_taurus.ARS-UCD1.2 \
#                 --fastqs=raw_data \
#                 --sample="$SAMPLE_ID"-1,"$SAMPLE_ID"-2,"$SAMPLE_ID"-3,"$SAMPLE_ID"-4 \
#                 --localcores=72 \
#                 --disable-ui

#mv "$SAMPLE_ID"/ ./cellranger_out/


cellranger count --id="$SAMPLE_ID" \
                 --transcriptome=Bos_taurus.ARS-UCD1.2 \
                 --fastqs=raw_data \
                 --sample="$SAMPLE_ID"-1,"$SAMPLE_ID"-2,"$SAMPLE_ID"-3,"$SAMPLE_ID"-4 \
                 --localcores=72 \
                 --disable-ui \
		 --force-cells=10000

mv "$SAMPLE_ID"/ ./cellranger_out_10000/
