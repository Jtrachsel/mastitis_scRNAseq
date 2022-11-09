#!/bin/bash

#SBATCH --job-name=fastqc
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 40
#SBATCH -t 48:00:00
#SBATCH -o ./logs/fastqc_"stdout.%j.%N"
#SBATCH -e ./logs/fastqc_"stderr.%j.%N"
#SBATCH --mem=64G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=

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
