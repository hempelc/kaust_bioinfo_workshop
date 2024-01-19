#### dada2 script from Larissa Fruehe ####
#### 09.03.2023 LF ####
#### Modified and polished by Chris Hempel
#### Used on IBEX supercomputer. When started with Rscript in the command line,
#### one parameter must be supplied at the end of the command, indicating the
#### number of CPUs to use. The script won't run if the number of CPUs is not
#### supplied. The number of CPUS for a SLURM run is stored in the variable
#### $SLURM_CPUS_PER_TASK
#### Example to use 32 CPUs: "Rscript --vanilla [name_of_this_script].R $SLURM_CPUS_PER_TASK"

# Check if number of CPUs has been supplied and stop script if not
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Supply number of CPUs as additional argument in the Rscript command", call. = FALSE)
}

# Load libraries
library(dada2)
library(ShortRead)
library(Biostrings)
library(ggplot2)
library(phyloseq)
library(plyr)

# Set parameters
## Full path to directory containing R1 and R2 sequences
sequencedir <- "/full/path/to/sequencedir/"
## Suffix of read files, for examples if reads are "XXX_R1/2-trimmed.fastq.gz", then suffix would be "-trimmed.fastq.gz"
## Note: names must contain R1/2
suffix <- "_trimmed.fastq.gz" # Replace as required
## Full path to output directory (will be created if non-existent)
outdir <- "/full/path/to/your/outdir"
## Filtering parameters for dada2
truncLen <- 140
maxEE <- 1
minLen <- 50
## Merging parameter for dada2
minOverlap <- 10
## Full paths to dada2 databases for taxonomy and species
taxonomy_db <- "/full/path/to/dada2_taxonomy.fasta"
species_db <- "/full/path/to/dada2_species.fasta"
## Specify if length filtering should be applied (TRUE if yes, FALSE if no)
length_filtering <- FALSE
## If length_filtering=TRUE, define upper and lower limits of expected amplicon length
amplicon_range_upper <- 294
amplicon_range_lower <- 212
## Number of CPUs (set automatically)
num_cpus <- as.numeric(args[1])


# Define a function to print datetime and text for log
print_datetime_and_text <- function(input_text) {
  # Get the current date and time
  current_date <- format(Sys.Date(), "%d.%m.%Y")
  current_time <- format(Sys.time(), "%H:%M")
  paste("===== ", current_date, " ", current_time, ": ", input_text, sep = "")
}
# Define a helper function to count the number of unique sequences
getN <- function(x) sum(getUniques(x))

print_datetime_and_text("Starting script.")

# Create output directories and change dir
dir.create(paste(outdir, "/filtered", sep = ""), recursive = TRUE, showWarnings = FALSE)
setwd(outdir)

# Retrieve a list of filenames of fastq.gz files and sort them
cutFs <- sort(list.files(sequencedir, pattern = paste0("R1", suffix), full.names = TRUE))
cutRs <- sort(list.files(sequencedir, pattern = paste0("R2", suffix), full.names = TRUE))

# Extract sample names from the filenames by splitting at "_R1" or "_R2", assumes filename = samplename_R1/R2.fastq.gz
sample.namesF <- sapply(strsplit(basename(cutFs), "_R1"), `[`, 1)
sample.namesR <- sapply(strsplit(basename(cutRs), "_R2"), `[`, 1)

# Check if the sample names from R1 and R2 files are identical
if (identical(sample.namesF, sample.namesR)) {
  print_datetime_and_text("Sanity check: forward and reverse read files are matching.....congratulations!")
} else {
  stop("Forward and reverse files do not match. Exiting script.")
}

# Assign names to the cutFs and cutRs vectors using the sample names
names(cutFs) <- sample.namesF
names(cutRs) <- sample.namesR

# Assign filenames to the filtered R1.fastq.gz and R2.fastq.gz files
filtFs <- file.path(paste(outdir, "/filtered", sep = ""), paste0(sample.namesF, "_filt_R1.fastq.gz"))
filtRs <- file.path(paste(outdir, "/filtered", sep = ""), paste0(sample.namesR, "_filt_R2.fastq.gz"))

# Filter and trim the cutFs and cutRs files using the specified parameters and store the result in out variable
print_datetime_and_text("Performing filtering and trimming...")
out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs,
  truncLen = truncLen, maxN = 0, maxEE = maxEE,
  truncQ = 2, minLen = minLen, rm.phix = TRUE,
  compress = TRUE, multithread = num_cpus
)
rownames(out) <- sapply(strsplit(basename(rownames(out)), "_R1"), `[`, 1)

# Save the out variable as Track.rds
saveRDS(out, "Track.rds")

# If some samples had zero reads after filtering, the dada fuction will give an error.
# Therefore, drop samples with 0 reads from the vectors
filtFs <- filtFs[file.exists(filtFs)]
filtRs <- filtRs[file.exists(filtRs)]
sample.namesF <- sapply(strsplit(basename(filtFs), "_filt"), `[`, 1)
sample.namesR <- sapply(strsplit(basename(filtRs), "_filt"), `[`, 1)

# Set a seed to ensure replicability of randomized steps
set.seed(100)

# Learn the error rates from the filtered fastq.gz files and save the result as errF/R.rds
print_datetime_and_text("Learning error rates...")
errF <- learnErrors(filtFs, multithread = num_cpus)
errR <- learnErrors(filtRs, multithread = num_cpus)
saveRDS(errF, "errF.rds")
saveRDS(errR, "errR.rds")

# Save a PDF file for the plot of forward and reverse error rates
pdf("Errorrates_fwd.pdf", width = 5, height = 5)
plotErrors(errF, nominalQ = TRUE)
dev.off()
pdf("Errorrates_rvs.pdf", width = 5, height = 5)
plotErrors(errR, nominalQ = TRUE)
dev.off()

# Uncomment the lines below if you want to read saved error rate objects
# errF <- readRDS("errF.rds")
# errR <- readRDS("errR.rds")

# Dereplicate reads
print_datetime_and_text("Dereplicating reads...")
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.namesF
names(derepRs) <- sample.namesR

# Run DADA2 algorithm on dereplicated reads
print_datetime_and_text("Running DADA2 algorithm on dereplicated reads...")
dadaFs <- dada(derepFs, err = errF, multithread = num_cpus, pool = "pseudo")
dadaRs <- dada(derepRs, err = errR, multithread = num_cpus, pool = "pseudo")
names(dadaFs) <- sample.namesF
names(dadaRs) <- sample.namesR
saveRDS(dadaFs, "dadaF.rds")
saveRDS(dadaRs, "dadaR.rds")

# Merge forward and reverse reads and save the result as mergers.rds
print_datetime_and_text("Merging forward and reverse reads...")
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, minOverlap = minOverlap, maxMismatch = 2, verbose = TRUE)
saveRDS(mergers, "mergers.rds")

# Create a sequence table from the merged pairs and save the result as seqtab.rds
seqtab <- makeSequenceTable(mergers)
saveRDS(seqtab, "seqtab.rds")

# Create a histogram of sequence lengths and save it as a JPEG file
## Define axis limits
axis_interval <- 5
low_len <- round_any(min(nchar(getSequences(seqtab))), 10, floor)
high_len <- round_any(max(nchar(getSequences(seqtab))), 10, ceiling)
range <- high_len - low_len
num_breaks <- range / axis_interval
## Make plot
jpeg(file = "Seqlength.jpg", res = 450, width = 15, height = 8, units = "in")
hist(nchar(getSequences(seqtab)), main = "Distribution of sequence lengths", breaks = range, xaxt = "n")
axis(1, at = seq(low_len, high_len, by = axis_interval))
dev.off()

# If length filtering activated, filter the seqtab based on the amplicon range and save the result
if (length_filtering) {
  print_datetime_and_text("Performing length filtering...")
  seqtab_filtered <- seqtab[, nchar(colnames(seqtab)) %in% seq(amplicon_range_lower, amplicon_range_upper)]
  saveRDS(seqtab_filtered, "seqtab_filtered.rds")

  # Remove chimeric sequences from the filtered sequence table and save the result
  print_datetime_and_text("Removing chimeras...")
  seqtab_nochim <- removeBimeraDenovo(seqtab_filtered, method = "consensus", multithread = num_cpus, verbose = TRUE)
  saveRDS(seqtab_nochim, "seqtab_nochim.rds")
} else {
  # Remove chimeric sequences from the sequence table and save the result
  print_datetime_and_text("Removing chimeras...")
  seqtab_nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = num_cpus, verbose = TRUE)
  saveRDS(seqtab_nochim, "seqtab_nochim.rds")
}

# Assign taxonomy to the non-chimeric sequences and save the result
print_datetime_and_text("Assigning taxonomy to the genus rank...")
taxa <- assignTaxonomy(seqtab_nochim, taxonomy_db, multithread = num_cpus)
print_datetime_and_text("Assigning taxonomy on the species rank...")
taxa <- addSpecies(taxa, species_db)
saveRDS(taxa, "Taxa.rds")

# Create a data frame with counts for various steps in the pipeline
if (length_filtering) {
  track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab_filtered), rowSums(seqtab_nochim))
  colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "seqtab", "seqtab_filt", "nonchim")
  track[track[, 2] == 0, 3:7] <- 0
} else {
  track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab_nochim))
  colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "seqtab", "nonchim")
  # Find rows where the filtered column contains 0 and replace values for following steps by 0
  track[track[, 2] == 0, 3:7] <- 0
}

# Write the data frame as "SeqOV.csv"
## Add header for rownames = samples
track$sample <- row.names(track)
row.names(track) <- NULL
## Reorder columns to have sample column as the first column
track <- track[, c("sample", names(track)[-which(names(track) == "sample")])]
write.csv(track, "SeqOV.csv")