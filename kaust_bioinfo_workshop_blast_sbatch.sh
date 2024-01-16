#!/bin/bash --login
#SBATCH --time=1:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=40
#SBATCH --partition=batch
#SBATCH -J kaust_bioinformatics_workshop_blastdbsetup
#SBATCH --output kaust_bioinformatics_workshop_blastdbsetup.%J.log

# Activate the blast module
module purge
module load blast

# Make the blast db
makeblastdb -in ~/kaust_bioinfo_workshop/kaust_bioinfo_workshop_files/blast_db_setup/ncbi_phosichthyidae_coi.fasta \
	-dbtype nucl
# Run blast
blastn -query ~/kaust_bioinfo_workshop/kaust_bioinfo_workshop_files/blast_db_setup/OTU_1.fasta \
 -db ~/kaust_bioinfo_workshop/kaust_bioinfo_workshop_files/blast_db_setup/ncbi_phosichthyidae_coi.fasta \
 -out blast_output.tsv \
 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxid" \
 -evalue 1e-05