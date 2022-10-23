#!/bin/bash

#SBATCH --job-name=reference_prep                              # name of the job submitted
#SBATCH -p short                                        # name of the queue you are submitting to
#SBATCH -N 1                                            # number of nodes in this job
#SBATCH -n 8                                           # number of cores/tasks in this job, you get all 20 cores with 2 threads per core with hyperthreading
#SBATCH -t 48:00:00                                     # time allocated for this job hours:mins:seconds
#SBATCH -o ./logs/reference_"stdout.%j.%N"                               # standard out %j adds job number to outputfile name and %N adds the node name
#SBATCH -e ./logs/reference_"stderr.%j.%N"                               # optional but it prints our standard error
#SBATCH --mem=64G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=julestrachsel@gmail.com



# modified from https://support.10xgenomics.com/single-cell-gene-expression/software/release-notes/build#hg19_3.0.0

# Genome metadata
genome="Bos_taurus.ARS-UCD1"
version="2"


# Set up source and build directories
build="Bos_taurus.ARS-UCD1.2_build"
mkdir -p "$build"


# Download source files if they do not exist in reference_sources/ folder
source="reference_sources"
mkdir -p "$source"


fasta_url="http://ftp.ensembl.org/pub/release-107/fasta/bos_taurus/dna/Bos_taurus.ARS-UCD1.2.dna.toplevel.fa.gz"
fasta_in="${source}/Bos_taurus.ARS-UCD1.2.dna.toplevel.fa"
gtf_url="http://ftp.ensembl.org/pub/release-107/gtf/bos_taurus/Bos_taurus.ARS-UCD1.2.107.gtf.gz"
gtf_in="${source}/Bos_taurus.ARS-UCD1.2.107.gtf"


if [ ! -f "$fasta_in" ]; then
    curl -sS "$fasta_url" | zcat > "$fasta_in"
fi
if [ ! -f "$gtf_in" ]; then
    curl -sS "$gtf_url" | zcat > "$gtf_in"
fi


cellranger mkgtf "$gtf_in" Bos_taurus.ARS-UCD1.2_build/Bos_taurus.ARS-UCD1.2.107_filtered.gtf \
	--attribute=gene_biotype:protein_coding \
	--attribute=gene_biotype:antisense \
	--attribute=gene_biotype:IG_C_gene \
	--attribute=gene_biotype:IG_C_pseudogene \
	--attribute=gene_biotype:IG_D_gene \
	--attribute=gene_biotype:IG_J_gene \
	--attribute=gene_biotype:IG_J_pseudogene \
	--attribute=gene_biotype:IG_LV_gene \
	--attribute=gene_biotype:IG_V_gene \
	--attribute=gene_biotype:IG_V_pseudogene \
	--attribute=gene_biotype:lincRNA \
	--attribute=gene_biotype:miRNA \
	--attribute=gene_biotype:misc_RNA \
	--attribute=gene_biotype:Mt_rRNA \
	--attribute=gene_biotype:pseudogene \
	--attribute=gene_biotype:ribozyme \
	--attribute=gene_biotype:rRNA \
	--attribute=gene_biotype:scaRNA \
	--attribute=gene_biotype:snoRNA \
	--attribute=gene_biotype:snRNA \
	--attribute=gene_biotype:sRNA \
	--attribute=gene_biotype:TR_C_gene \
	--attribute=gene_biotype:TR_D_gene \
	--attribute=gene_biotype:TR_J_gene \
	--attribute=gene_biotype:TR_J_pseudogene \
	--attribute=gene_biotype:TR_V_gene \
	--attribute=gene_biotype:TR_V_pseudogene

cellranger mkref --genome=Bos_taurus.ARS-UCD1.2 \
                 --fasta=reference_sources/Bos_taurus.ARS-UCD1.2.dna.toplevel.fa \
                 --genes=Bos_taurus.ARS-UCD1.2_build/Bos_taurus.ARS-UCD1.2.107_filtered.gtf \
                 --ref-version=1.2.107 \
		 --nthreads=8



#sbatch scripts/03_fastQC.slurm

