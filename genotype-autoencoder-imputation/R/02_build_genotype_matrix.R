# ------------------------------------------------------------------
# 02_build_genotype_matrix.R
#
# Stage 2: Build the full SNP x sample dosage matrix from IMPUTE2
# output. This matrix is the input the autoencoder trains on.
#
# Output: outputs/toy_genotype_matrix.csv
#   rows    = samples
#   columns = SNPs (dosage values, 0-2)
# ------------------------------------------------------------------

library(data.table)
library(yaml)

cfg <- read_yaml("config.yml")
base_directory <- cfg$paths$base_directory

imputed <- fread(cfg$paths$output_prefix, header = FALSE)
# imputed columns: ID, rsID, pos, A1, A2, then 3 probabilities per sample

sample_file <- file.path(base_directory, cfg$paths$gene_files$sample_file)
sample_df <- fread(sample_file, header = TRUE)
sample_ids <- sample_df$ID_2
if (sample_ids[1] == "0") sample_ids <- sample_ids[-1]

n_samples <- length(sample_ids)

get_dosage_from_row <- function(row, n_samples) {
  row_vec <- as.numeric(row)
  probs   <- row_vec[6:(5 + 3 * n_samples)]
  prob_mat <- matrix(probs, ncol = 3, byrow = TRUE)
  as.numeric(prob_mat %*% c(0, 1, 2))
}

# Build dosage matrix across all SNPs (not just TAS2R38)
dos_mat <- t(apply(imputed, 1, get_dosage_from_row, n_samples = n_samples))
colnames(dos_mat) <- sample_ids   # columns = samples
rownames(dos_mat) <- imputed$V1   # rows = SNP IDs

# Transpose so rows = samples, columns = SNPs (autoencoder expects this shape)
X <- t(dos_mat)

write.csv(X, file = cfg$paths$genotype_matrix_csv, row.names = TRUE)
