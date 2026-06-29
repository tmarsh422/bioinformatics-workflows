# ------------------------------------------------------------------
# 03_compare_results.R
#
# Stage 4: Compare autoencoder-imputed dosages against IMPUTE2's
# dosages for the TAS2R38 SNPs. Computes correlation, RMSE, and
# hard-call genotype accuracy, and produces diagnostic plots.
# ------------------------------------------------------------------

library(data.table)
library(ggplot2)
library(yaml)

cfg <- read_yaml("config.yml")
tas_snps <- cfg$target$snps

# ---- Load both sets of dosages ------------------------------------

tas_imp <- fread(cfg$paths$impute2_dosages_csv)
ae_raw  <- fread(cfg$paths$ae_dosages_csv)

ae_tas <- ae_raw[, c("sample_id", tas_snps), with = FALSE]

# ---- Align samples between the two sources -------------------------

common_ids <- intersect(tas_imp$sample_id, ae_tas$sample_id)

tas_imp_f <- tas_imp[sample_id %in% common_ids]
ae_tas_f  <- ae_tas[sample_id %in% common_ids]

setkey(tas_imp_f, sample_id)
setkey(ae_tas_f,  sample_id)
tas_imp_f <- tas_imp_f[order(sample_id)]
ae_tas_f  <- ae_tas_f[order(sample_id)]

stopifnot(identical(tas_imp_f$sample_id, ae_tas_f$sample_id))

# ---- Reshape to long format for comparison -------------------------

imp_long <- melt(tas_imp_f, id.vars = "sample_id", variable.name = "snp", value.name = "dos_imp")
ae_long  <- melt(ae_tas_f,  id.vars = "sample_id", variable.name = "snp", value.name = "dos_ae")

comp <- merge(imp_long, ae_long, by = c("sample_id", "snp"))
comp[, diff := dos_ae - dos_imp]

# ---- Continuous metrics: correlation & RMSE ------------------------

rmse_overall <- sqrt(mean(comp$diff^2, na.rm = TRUE))
cor_overall  <- cor(comp$dos_imp, comp$dos_ae, use = "complete.obs")

cat("Overall RMSE:", rmse_overall, "\n")
cat("Overall correlation:", cor_overall, "\n")

comp_stats <- comp[, .(
  cor  = cor(dos_imp, dos_ae, use = "complete.obs"),
  rmse = sqrt(mean((dos_ae - dos_imp)^2, na.rm = TRUE))
), by = snp]

print(comp_stats)

# ---- Hard-call genotype accuracy (0/1/2) ----------------------------

to_geno <- function(d) {
  ifelse(d < 0.5, 0L, ifelse(d < 1.5, 1L, 2L))
}

comp[, geno_imp := to_geno(dos_imp)]
comp[, geno_ae  := to_geno(dos_ae)]

acc_overall <- mean(comp$geno_imp == comp$geno_ae, na.rm = TRUE)
cat("Overall genotype accuracy:", acc_overall, "\n")

acc_by_snp <- comp[, .(accuracy = mean(geno_imp == geno_ae, na.rm = TRUE)), by = snp]
print(acc_by_snp)

# ---- Plots -----------------------------------------------------------

target_snp_focus <- cfg$target$hide_snp

p_scatter <- ggplot(comp[snp == target_snp_focus], aes(x = dos_imp, y = dos_ae)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0) +
  labs(
    title = paste("IMPUTE2 vs Autoencoder Dosages for", target_snp_focus),
    x = "IMPUTE2 dosage",
    y = "Autoencoder dosage"
  )
print(p_scatter)

p_bland_altman <- ggplot(comp, aes(x = (dos_imp + dos_ae) / 2, y = diff, color = snp)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(alpha = 0.5) +
  labs(
    title = "Difference: Autoencoder - IMPUTE2 Dosages",
    x = "Mean dosage (IMPUTE2 & AE)",
    y = "Difference (AE - IMPUTE2)"
  )
print(p_bland_altman)

sub <- comp[snp %in% tas_snps]
sub[, err_ae := abs(dos_imp - dos_ae)]

p_error <- ggplot(sub, aes(x = snp, y = err_ae)) +
  geom_point(alpha = 0.6) +
  labs(
    title = paste("Imputation Error for", target_snp_focus),
    subtitle = "Autoencoder vs IMPUTE2",
    x = "SNP identity",
    y = "Absolute error (autoencoder)"
  ) +
  theme_minimal()
print(p_error)
