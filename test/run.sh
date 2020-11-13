#!/bin/bash
set -e

module load singularity
eval "$(conda shell.bash hook)"
conda activate snakemake

cd ..
snakemake --cores 1 --configfile test/config/test.yml --use-singularity
