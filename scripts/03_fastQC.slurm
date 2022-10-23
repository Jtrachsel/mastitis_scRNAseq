#!/bin/bash

#SBATCH --job-name=fastqc                               # name of the job submitted
#SBATCH -p short                                        # name of the queue you are submitting to
#SBATCH -N 1                                            # number of nodes in this job
#SBATCH -n 40                                           # number of cores/tasks in this job, you get all 20 cores with 2 threads per core with hyperthreading
#SBATCH -t 48:00:00                                     # time allocated for this job hours:mins:seconds
#SBATCH -o ./logs/fastqc_"stdout.%j.%N"                               # standard out %j adds job number to outputfile name and %N adds the node name
#SBATCH -e ./logs/fastqc_"stderr.%j.%N"                               # optional but it prints our standard error
#SBATCH --mem=64G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=julestrachsel@gmail.com

module load fastqc
module load miniconda
source activate /project/fsepru/conda_envs/multiqc


if [ ! -d fastqc_results/ ]; then
	
	mkdir -p fastqc_results
	mkdir -p fastqc_results/lane4
	mkdir -p fastqc_results/lane5
	mkdir -p fastqc_results/lane6

	fastqc -o fastqc_results/lane4 -t 40 raw_data/lane4/*fastq
	fastqc -o fastqc_results/lane5 -t 40 raw_data/lane5/*fastq
	fastqc -o fastqc_results/lane6 -t 40 raw_data/lane6/*fastq

	multiqc -n multiQC_lane4.html fastqc_results/lane4
	multiqc -n multiQC_lane5.html fastqc_results/lane5
	multiqc -n multiQC_lane6.html fastqc_results/lane6

	mv multiQC_lane*/ ./fastqc_results
	mv multiQC*html ./fastqc_results
fi


# everything looks good. maybe one sample that has issues.

sbatch scripts/04_cellranger_map_individual.slurm
