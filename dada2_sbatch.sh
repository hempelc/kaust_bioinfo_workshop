#!/bin/bash
#SBATCH --job-name=dada2_ADD_PROJECT_NAME
#SBATCH --output=dada2_ADD_PROJECT_NAME.%J.log
#SBATCH --partition=batch
#SBATCH --cpus-per-task=40
#SBATCH --time=48:00:00 # Adjust as necessary via trial and error
#SBATCH --mem=64G # Adjust as necessary via trial and error

module load R

# Start the R script and specify the number of CPUs used (it's automatically stored in the variable $SLURM_CPUS_PER_TASK, because we used the cpus-per-task SBATCH command)
Rscript --vanilla dada2_master_script.R $SLURM_CPUS_PER_TASK # If you changed the name of the R script, it needs to be changed here too
