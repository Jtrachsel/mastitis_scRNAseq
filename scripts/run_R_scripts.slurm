#!/bin/bash

#SBATCH --job-name=scRNA_pipe
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -t 48:00:00
#SBATCH -o "stdout.%j.%N"
#SBATCH -e "stderr.%j.%N"
#SBATCH --mem=200G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=julestrachsel@gmail.com

set -e
module load singularity

## need to make sure renv is set up
# https://stackoverflow.com/questions/66311017/rscript-to-use-renv-environment

# Ensure your working directory is set to the root of your renv project,
# and that the renv project's auto-loader is active.
# You can set up the auto-loader by calling renv::activate() from R in that project.

# I did this prior to running this script.  Maybe I should include a 'setup' script.

#
# I developed these scripts using RStudio via open on demand on Ceres
# this uses a singularity container for the RStudio server instance
# this means that all the packages I installed during this process
# depend on the singularity container libraries.
# So we have to use that same container to run the R scripts in a non-interactive way

#singularity run /reference/containers/RStudio/4.2.0/rocker-geospatial.sif Rscript scripts/05_check_mapping.R
#singularity run /reference/containers/RStudio/4.2.0/rocker-geospatial.sif Rscript scripts/06_map_gene_ids.R
#singularity run /reference/containers/RStudio/4.2.0/rocker-geospatial.sif Rscript scripts/07_SoupX.R
#singularity run /reference/containers/RStudio/4.2.0/rocker-geospatial.sif Rscript scripts/08_scDblFinder.R

# this takes ~ 1 min
#date
#echo 09_QC.R
#singularity run /reference/containers/RStudio/4.2.0/rocker-geospatial.sif Rscript scripts/09_QC.R

#date
#echo 10_integrate.R
#singularity run /reference/containers/RStudio/4.2.0/rocker-geospatial.sif Rscript scripts/10_integrate.R

#date
#echo 11_cluster.R
#singularity run /reference/containers/RStudio/4.2.0/rocker-geospatial.sif Rscript scripts/11_cluster.R

date
echo 12_marker_genes.R
singularity run /reference/containers/RStudio/4.2.0/rocker-geospatial.sif Rscript scripts/12_marker_genes.R

date
echo 13_general_cluster_identity.R
singularity run /reference/containers/RStudio/4.2.0/rocker-geospatial.sif Rscript scripts/13_general_cluster_identity.R
