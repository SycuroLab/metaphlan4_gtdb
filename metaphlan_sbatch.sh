#!/bin/bash

#SBATCH --partition=synergy
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=2-00:00:00
#SBATCH --mem=5G
#SBATCH --error=run_metaphlan_sbatch.%J.err
#SBATCH --output=run_metaphlan_sbatch.%J.out

log_dir="$(pwd)"
log_file="logs/metaphlan-analysis.log.txt"
num_jobs=200

echo "started at: `date`"

# Load the ~/.bashrc file as source.
source ~/.bashrc

# Activate the snakemake conda environment.
conda activate snakemake

# Unlock snakemake folder as a fail safe.
snakemake --unlock

snakemake --latency-wait 25 --rerun-incomplete --cluster-config cluster.json --cluster 'sbatch --partition={cluster.partition} --cpus-per-task={cluster.cpus-per-task} --nodes={cluster.nodes} --ntasks={cluster.ntasks} --time={cluster.time} --mem={cluster.mem} --output={cluster.output} --error={cluster.error}' --jobs $num_jobs --use-conda &> $log_dir/$log_file

echo "finished with exit code $? at: `date`"


