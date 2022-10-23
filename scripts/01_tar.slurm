#!/bin/bash

#SBATCH --job-name=tar                              # name of the job submitted
#SBATCH -p short                                        # name of the queue you are submitting to
#SBATCH -N 1                                            # number of nodes in this job
#SBATCH -n 40                                           # number of cores/tasks in this job, you get all 20 cores with 2 threads per core with hyperthreading
#SBATCH -t 48:00:00                                     # time allocated for this job hours:mins:seconds
#SBATCH -o "stdout.%j.%N"                               # standard out %j adds job number to outputfile name and %N adds the node name
#SBATCH -e "stderr.%j.%N"                               # optional but it prints our standard error
#SBATCH --mem=64G   
#SBATCH --mail-type=ALL
#SBATCH --mail-user=julestrachsel@gmail.com

set -e
module load parallel

# ENTER COMMANDS HERE:

# stage project on 90daydata partition

mkdir -p logs
mkdir -p /90daydata/fsep_004/john_lippolis_milk_blood/
cp -r -n * /90daydata/fsep_004/john_lippolis_milk_blood/

cd /90daydata/fsep_004/john_lippolis_milk_blood/

# untar the archive

if [ ! -d ./raw_data/lane4/ ]; then

	tar -xvf "./raw_data/4,5,6_31025_01_1312milk_H7LMV_1705.tar" -C ./raw_data/

	parallel 'gunzip {}' ::: ./raw_data/lane4/*gz
	parallel 'gunzip {}' ::: ./raw_data/lane5/*gz
	parallel 'gunzip {}' ::: ./raw_data/lane6/*gz

fi

# generate a list of sample names from the fastqs
mkdir -p outputs
for x in raw_data/lane4/*R1*; do y=$(basename $x); echo "${y%_S*}"; done | grep -v Undetermined > outputs/sample_names.txt

# extract some metadata from the sample names
sed -E 's/([0-9]+)([a-z]+)-([0-9]+)/&\t\1\t\2\t\3/g' outputs/sample_names.txt > outputs/metadata.tsv

# get biological sample names from replicates
cat outputs/sample_names.txt |sed 's/-[0-9]//g' |uniq > outputs/bio_samples.txt

# start the next script


sbatch scripts/02_cellranger_reference_prep_script.slurm


#End of file
