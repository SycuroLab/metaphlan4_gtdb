#!/bin/bash

#SBATCH --partition=synergy,cpu2019,cpu2021,cpu2022,cpu2023
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=14
#SBATCH --time=7-00:00:00
#SBATCH --mem=50G
#SBATCH --error=run_install_metaphlan4_database.%J.err
#SBATCH --output=run_install_metaphlan4_database.%J.out

# Load the ~/.bashrc file as source.
source ~/.bashrc

# Activate the metaphlan4 conda environment.
conda activate metaphlan4_env

# The metaphlan database directory.
metaphlan_db_dir="/bulk/IMCshared_bulk/shared/dbs/metaphlan4.0.5"

# The number of threads to use.
num_threads=14

# Install the metaphlan4 database.
metaphlan --install --nproc ${num_threads} --bowtie2db ${metaphlan_db_dir}

