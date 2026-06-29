# ------------------------------------------------------------------
# 01_run_impute2.R
#
# Stage 1: Run IMPUTE2 on chr7 genotype data, extract imputed dosages
# for the three TAS2R38 SNPs, and write them to a CSV for later
# comparison against the autoencoder's predictions.
#
# Reads all paths/parameters from config.yml — no hardcoded paths.
# ------------------------------------------------------------------

library(data.table)
library(tidyverse)
library(yaml)

cfg <- read_yaml("config.yml")
base_directory <- cfg$paths$base_directory

impute2_exe <- file.path(base_directory, cfg$paths$impute2_exe)
map_file    <- file.path(base_directory, cfg$paths$ref_panels$map_file)
hap_file    <- file.path(base_directory, cfg$paths$ref_panels$hap_file)
legend_file <- file.path(base_directory, cfg$paths$ref_panels$legend_file)

gens_file   <- file.path(base_directory, cfg$paths$gene_files$gens_file)
sample_file <- file.path(base_directory, cfg$paths$gene_files$sample_file)
strand_file <- file.path(base_directory, cfg$paths$gene_files$strand_file)

output_prefix <- cfg$paths$output_prefix
tas_snps <- cfg$target$snps

# ---- Run IMPUTE2 --------------------------------------------------

impute_args <- c(
  "-m", shQuote(map_file),               # SNP positional info
  "-h", shQuote(hap_file),               # reference haplotypes
  "-l", shQuote(legend_file),            # variants in study sample
  "-g", shQuote(gens_file),              # genotype probabilities per SNP
  "-sample_g", shQuote(sample_file),     # study sample list
  "-strand_g", shQuote(strand_file),
  "-Ne", as.character(cfg$impute2$Ne),
  "-int", as.character(cfg$impute2$interval[1]), as.character(cfg$impute2$interval[2]),
  "-o", shQuote(output_prefix)
)

system2(impute2_exe, args = impute_args)

# ---- Inspect imputation quality for target SNPs -------------------

info_file <- paste0(output_prefix, "_info")
info <- read.table(info_file, header = TRUE, stringsAsFactors = FALSE)
info_tas <- subset(info, snp_id %in% tas_snps)
print(info_tas)

# ---- Extract dosages for the target SNPs --------------------------
# Dosage = expected allele count at a locus (continuous value, 0-2)

imputed <- fread(output_prefix, header = FALSE)
imputed_tas <- imputed[V1 %in% tas_snps]

n_cols <- ncol(imputed_tas)
n_samples <- (n_cols - 5) / 3

sample_df <- read.table(sample_file, header = TRUE, stringsAsFactors = FALSE)
sample_ids <- sample_df$ID_2
if (sample_ids[1] == "0") sample_ids <- sample_ids[-1]

get_dosage_from_row <- function(row, n_samples) {
  row_vec <- as.numeric(row)
  probs <- row_vec[6:(5 + 3 * n_samples)]
  prob_mat <- matrix(probs, ncol = 3, byrow = TRUE)
  # dosage = 0 * P(AA) + 1 * P(AB) + 2 * P(BB)
  as.numeric(prob_mat %*% c(0, 1, 2))
}

dosage_list <- lapply(tas_snps, function(snp) {
  get_dosage_from_row(imputed_tas[V1 == snp], n_samples)
})

tas_dosages <- data.frame(sample_id = sample_ids)
for (i in seq_along(tas_snps)) {
  tas_dosages[[tas_snps[i]]] <- dosage_list[[i]]
}

head(tas_dosages)

write.csv(tas_dosages, cfg$paths$impute2_dosages_csv, row.names = FALSE)
