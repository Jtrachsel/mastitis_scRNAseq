!/bin/bash

#SBATCH --job-name=soupx                             
#SBATCH -p short                                 
#SBATCH -N 1                                     
#SBATCH -n 8                                    
#SBATCH -t 48:00:00                              
#SBATCH -o "stdout.%j.%N"                        
#SBATCH -e "stderr.%j.%N"                        
#SBATCH --mem=64G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=julestrachsel@gmail.com

set -e
module load singularity

singularity run /reference/containers/RStudio/4.2.0/rocker-geospatial.sif Rscript ./scripts/soupx.R

