#!/bin/bash --login
#SBATCH --time=3:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=40
#SBATCH --partition=batch
#SBATCH -J kaust_bioinformatics_workshop_apscale
#SBATCH --output kaust_bioinformatics_workshop_apscale.%J.log

# Activate the apscale mamba environment
module purge
mamba deactivate
mamba activate apscale

# Load required module blast
module load blast

# Run apscale
apscale_wrapper.py --sequence_dir example_sequences \
	--project_name kaust_bioinformatics_workshop \
	--forward_primer GGWACWGGWTGAACWGTWTAYCCYCC \
	--reverse_primer TANACYTCNGGRTGNCCRAARAAYCA \
	--min_length 303 \
	--max_length 323 \
	--clusteringtool swarm \
	--coi True \
	--prior_denoising True \
	--remove_negative_controls True \
	--negative_controls LF-STNEG1 \
	--run_blast True \
	--database MIDORI2_UNIQ_NUC_GB259_CO1_BLAST/MIDORI2_UNIQ_NUC_GB259_CO1_BLAST \
	--database_format midori2 \
	--cores 40
