#!/bin/bash

#SBATCH --job-name='saige'
#SBATCH --output=saige-%j-stdout.log
#SBATCH --error=saige-%j-stderr.log
#SBATCH --cpus-per-task=9
#SBATCH --mem=300G
#SBATCH --partition=long

#Run command

nextflow run main.nf -c conf/test.config -profile singularity -resume
