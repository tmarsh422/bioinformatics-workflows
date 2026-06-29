[README (1).md](https://github.com/user-attachments/files/29477846/README.1.md)
# Genotype Imputation: Autoencoder vs. IMPUTE2

A small research pipeline comparing a PyTorch denoising autoencoder
against [IMPUTE2](https://mathgen.stats.ox.ac.uk/impute/impute_v2.html)
— a standard statistical-genetics imputation tool — for predicting
genotype dosages at three SNPs in the *TAS2R38* bitter-taste receptor
gene.

## Why

IMPUTE2 imputes missing genotypes using a reference haplotype panel
(1000 Genomes) and a population-genetics model of linkage
disequilibrium (LD). This project asks: **can a simple neural network
learn comparable LD structure directly from the study sample's own
genotypes, without a reference panel?** A denoising autoencoder is
trained to reconstruct a held-out SNP from the dosages of all other
genotyped SNPs in the window, and its predictions are compared against
IMPUTE2's on correlation, RMSE, and hard-call genotype accuracy.

## Pipeline

```
IMPUTE2 (R) ──> dosage extraction (R) ──> genotype matrix (R)
                                                │
                                                ▼
                                  autoencoder training (Python)
                                                │
                                                ▼
                                  comparison & plots (R)
```

| Stage | Script | What it does |
|---|---|---|
| 1 | `R/01_run_impute2.R` | Runs IMPUTE2, extracts TAS2R38 dosages |
| 2 | `R/02_build_genotype_matrix.R` | Builds the full sample × SNP dosage matrix |
| 3 | `python/train.py` (+ `python/autoencoder.py`) | Trains the autoencoder, imputes the target SNP |
| 4 | `R/03_compare_results.R` | Computes correlation/RMSE/accuracy, generates plots |

A narrative walkthrough combining all four stages with discussion is
in [`notebooks/autoencoder_walkthrough.qmd`](notebooks/autoencoder_walkthrough.qmd).

## Project structure

```
.
├── config.yml                          # all paths & hyperparameters
├── environment.yml                     # conda environment (R + Python)
├── data/
│   └── README.md                       # expected data layout (not committed)
├── R/
│   ├── 01_run_impute2.R
│   ├── 02_build_genotype_matrix.R
│   └── 03_compare_results.R
├── python/
│   ├── autoencoder.py                  # model + train/inference logic
│   └── train.py                        # CLI entry point
├── notebooks/
│   └── autoencoder_walkthrough.qmd     # narrative version with discussion
└── outputs/                            # generated CSVs & figures (not committed)
```

## Setup

```bash
conda env create -f environment.yml
conda activate genotype-ae
```

Then populate `data/raw/` per [`data/README.md`](data/README.md) and
adjust `config.yml` if your file layout differs.

## Running

```bash
# Stage 1-2: IMPUTE2 + genotype matrix construction
Rscript R/01_run_impute2.R
Rscript R/02_build_genotype_matrix.R

# Stage 3: train the autoencoder and impute the target SNP
python python/train.py

# Stage 4: compare results
Rscript R/03_compare_results.R
```

Or render the full walkthrough:

```bash
quarto render notebooks/autoencoder_walkthrough.qmd
```

## Method notes

- **Dosage**: the expected allele count at a locus, a continuous value
  in `[0, 2]` rather than a hard genotype call.
- **Training signal**: in each batch, the target SNP is randomly
  masked for ~30% of samples (`p_hide` in `config.yml`), and the
  reconstruction loss is computed only at those masked positions —
  this is what forces the network to predict the SNP from the LD
  structure of the others rather than trivially copying its input.
- **Evaluation**: correlation, RMSE, and 0/1/2 hard-call accuracy
  against IMPUTE2's dosages, plus a Bland-Altman plot to check for
  systematic bias.

## Results

_Add a short summary here once you've run the pipeline on your data —
e.g. overall correlation/RMSE, and whether the autoencoder over- or
under-performs IMPUTE2 for each SNP. A figure or two from
`outputs/` makes this section land well for anyone skimming the repo._

## Limitations

The autoencoder has no access to an external reference panel — it can
only exploit LD signal present in the genotyped sample itself, so its
agreement with IMPUTE2 is a rough concordance check, not a
ground-truth validation. Performance is also expected to be sensitive
to sample size and the density of SNPs genotyped in the surrounding
window.
