#!/bin/bash

#SBATCH --job-name=reference_prep                              # name of the job submitted
#SBATCH -p short                                        # name of the queue you are submitting to
#SBATCH -N 1                                            # number of nodes in this job
#SBATCH -n 4                                           # number of cores/tasks in this job, you get all 20 cores with 2 threads per core with hyperthreading
#SBATCH -t 48:00:00                                     # time allocated for this job hours:mins:seconds
#SBATCH -o "stdout.%j.%N"                               # standard out %j adds job number to outputfile name and %N adds the node name
#SBATCH -e "stderr.%j.%N"                               # optional but it prints our standard error
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

### Jules Note
# did not modify the fasta

# Modify sequence headers in the Ensembl FASTA to match the file
# "GRCh38.primary_assembly.genome.fa" from GENCODE. Unplaced and unlocalized
# sequences such as "KI270728.1" have the same names in both versions.
#
# Input FASTA:
#   >1 dna:chromosome chromosome:GRCh38:1:1:248956422:1 REF
#
# Output FASTA:
#   >chr1 1
fasta_modified="$build/$(basename "$fasta_in").modified"
# sed commands:
# 1. Replace metadata after space with original contig name, as in GENCODE
# 2. Add "chr" to names of autosomes and sex chromosomes
# 3. Handle the mitochrondrial chromosome
#cat "$fasta_in" \
#    | sed -E 's/^>(\S+).*/>\1 \1/' \
#    | sed -E 's/^>([0-9]+|[XY]) />chr\1 /' \
#    | sed -E 's/^>MT />chrM /' \
#    > "$fasta_modified"



### Jules Note: 
###   It seems like the Bos taurus genome version we are using doesn't need this step
###   no suffixes to remove 

# Remove version suffix from transcript, gene, and exon IDs in order to match
# previous Cell Ranger reference packages
#
# Input GTF:
#     ... gene_id "ENSG00000223972.5"; ...
# Output GTF:
#     ... gene_id "ENSG00000223972"; gene_version "5"; ...
#gtf_modified="$build/$(basename "$gtf_in").modified"
# Pattern matches Ensembl gene, transcript, and exon IDs for human or mouse:
#ID="(ENS(BTA)?[GTE][0-9]+)\.([0-9]+)"
#cat "$gtf_in" \
#    | sed -E 's/gene_id "'"$ID"'";/gene_id "\1"; gene_version "\3";/' \
#    | sed -E 's/transcript_id "'"$ID"'";/transcript_id "\1"; transcript_version "\3";/' \
#    | sed -E 's/exon_id "'"$ID"'";/exon_id "\1"; exon_version "\3";/' \
#    > "$gtf_modified"




# Define string patterns for GTF tags
# NOTES:
# - Since GENCODE release 31/M22 (Ensembl 97), the "lincRNA" and "antisense"
#   biotypes are part of a more generic "lncRNA" biotype.
# - These filters are relevant only to GTF files from GENCODE. The GTFs from
#   Ensembl release 98 have the following differences:
#   - The names "gene_biotype" and "transcript_biotype" are used instead of
#     "gene_type" and "transcript_type".
#   - Readthrough transcripts are present but are not marked with the
#     "readthrough_transcript" tag.
#   - Only the X chromosome versions of genes in the pseudoautosomal regions
#     are present, so there is no "PAR" tag.

### Jules Note
# I had to change 'gene_type' to 'gene_biotype'
# and 'transcript_type' to 'transcript_biotype' below

BIOTYPE_PATTERN=\
"(protein_coding|lncRNA|\
IG_C_gene|IG_D_gene|IG_J_gene|IG_LV_gene|IG_V_gene|\
IG_V_pseudogene|IG_J_pseudogene|IG_C_pseudogene|\
TR_C_gene|TR_D_gene|TR_J_gene|TR_V_gene|\
TR_V_pseudogene|TR_J_pseudogene)"
GENE_PATTERN="gene_biotype \"${BIOTYPE_PATTERN}\""
TX_PATTERN="transcript_biotype \"${BIOTYPE_PATTERN}\""
READTHROUGH_PATTERN="tag \"readthrough_transcript\"" # NO OCCURANCE OF 'readthrough' 
PAR_PATTERN="tag \"PAR\""                            # no occurange of 'PAR' in tags


# Construct the gene ID allowlist. We filter the list of all transcripts
# based on these criteria:
#   - allowable gene_type (biotype)
#   - allowable transcript_type (biotype)
#   - no "PAR" tag (only present for Y chromosome PAR)
#   - no "readthrough_transcript" tag
# We then collect the list of gene IDs that have at least one associated
# transcript passing the filters.

### Jules note
# since we commended out the first modification to the gtf file, I have changed
# the "$gtf_modified" to "$gtf_in" here and below

# cat "$gtf_modified" \
cat "$gtf_in" \
    | awk '$3 == "transcript"' \
    | grep -E "$GENE_PATTERN" \
    | grep -E "$TX_PATTERN" \
    | grep -Ev "$READTHROUGH_PATTERN" \
    | grep -Ev "$PAR_PATTERN" \
    | sed -E 's/.*(gene_id "[^"]+").*/\1/' \
    | sort \
    | uniq \
    > "${build}/gene_allowlist"

# how many?


# Filter the GTF file based on the gene allowlist
gtf_filtered="${build}/$(basename "$gtf_in").filtered"
# Copy header lines beginning with "#"
grep -E "^#" "$gtf_in" > "$gtf_filtered"
# Filter to the gene allowlist
grep -Ff "${build}/gene_allowlist" "$gtf_in" \
    >> "$gtf_filtered"

# final number of entries?

# Create reference package
cellranger mkref --ref-version="$version" \
    --genome="$genome" --fasta="$fasta_in" --genes="$gtf_filtered"


# start the next script
sbatch scripts/03_fastQC.slurm

