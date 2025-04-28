#!/bin/bash
#SBATCH --job-name=st-pipeline
#SBATCH --output=logs_slurm/st-pipeline_%j.out
#SBATCH --error=logs_slurm/st-pipeline_%j.err
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=250G

# Create logs directory if it doesn't exist
mkdir -p logs_slurm

# Activate st conda environment
source ~/.bashrc
conda activate st

# Run snakemake rule all
snakemake --cores all all_annotation \
    --use-conda \
    --rerun-incomplete \
    --keep-going \
    --latency-wait 60 \
    --keep-incomplete \
#    --ignore-incomplete

# Deactivate conda environment
conda deactivate 


