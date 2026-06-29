"""
train.py

Entry point: loads the genotype matrix produced by the R pipeline,
trains the autoencoder to impute a target SNP, and writes the
resulting dosages to CSV for comparison against IMPUTE2.

Usage:
    python python/train.py
"""

import yaml
import numpy as np
import pandas as pd
import torch

from autoencoder import (
    AEConfig,
    normalize,
    build_missing_mask,
    train_autoencoder,
    impute_target_snp,
)


def load_config(path: str = "config.yml") -> dict:
    with open(path) as f:
        return yaml.safe_load(f)


def main():
    cfg = load_config()

    df = pd.read_csv(cfg["paths"]["genotype_matrix_csv"], index_col=0)
    sample_ids = df.index.to_list()
    snp_ids = df.columns.to_list()

    X_np = df.to_numpy().astype(np.float32)

    X_norm, snp_mean, snp_std = normalize(X_np)
    missing_mask = build_missing_mask(X_norm)
    X_norm_filled = np.nan_to_num(X_norm, nan=0.0)

    target_snp = cfg["target"]["hide_snp"]
    target_j = snp_ids.index(target_snp)

    ae_cfg = AEConfig(**cfg["autoencoder"])

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    model = train_autoencoder(
        X_norm_filled, missing_mask, target_j, ae_cfg, device
    )

    recon_dosages = impute_target_snp(
        model, X_norm_filled, target_j, snp_mean, snp_std, device
    )

    ae_full_df = pd.DataFrame(recon_dosages, index=sample_ids, columns=snp_ids)

    tas_snps = cfg["target"]["snps"]
    ae_dos_df = ae_full_df[tas_snps].copy()
    ae_dos_df.insert(0, "sample_id", sample_ids)

    ae_dos_df.to_csv(cfg["paths"]["ae_dosages_csv"], index=False)
    print(f"Wrote autoencoder dosages to {cfg['paths']['ae_dosages_csv']}")


if __name__ == "__main__":
    main()
