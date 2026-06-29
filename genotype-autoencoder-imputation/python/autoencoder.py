"""
autoencoder.py

Denoising autoencoder for SNP genotype dosage imputation.

The model is trained on a full sample x SNP dosage matrix. During
training, a target SNP is randomly hidden for a fraction of samples
in each batch, and the loss is computed only on the hidden positions
-- forcing the model to learn to reconstruct that SNP from the
genotype context of the others (linkage disequilibrium structure).
"""

from dataclasses import dataclass

import numpy as np
import torch
from torch import nn
from torch.utils.data import DataLoader, TensorDataset


@dataclass
class AEConfig:
    bottleneck: int = 64
    hidden_dim: int = 256
    p_hide: float = 0.30
    batch_size: int = 64
    n_epochs: int = 265
    learning_rate: float = 1e-3


class AE(nn.Module):
    """Simple fully-connected autoencoder."""

    def __init__(self, n_features: int, bottleneck: int = 64, hidden_dim: int = 256):
        super().__init__()
        self.encoder = nn.Sequential(
            nn.Linear(n_features, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, bottleneck),
            nn.ReLU(),
        )
        self.decoder = nn.Sequential(
            nn.Linear(bottleneck, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, n_features),
        )

    def forward(self, x):
        z = self.encoder(x)
        return self.decoder(z)


def normalize(X: np.ndarray):
    """Per-SNP z-score normalization. Returns (X_norm, mean, std)."""
    mean = X.mean(axis=0, keepdims=True)
    std = X.std(axis=0, keepdims=True) + 1e-6
    return (X - mean) / std, mean, std


def build_missing_mask(X_norm: np.ndarray):
    """True where a value is observed (not NaN)."""
    return ~np.isnan(X_norm)


def train_autoencoder(
    X_norm_filled: np.ndarray,
    missing_mask: np.ndarray,
    target_j: int,
    cfg: AEConfig,
    device: torch.device,
    verbose: bool = True,
) -> AE:
    """
    Train the autoencoder to impute the SNP at column index `target_j`.

    For each batch, a random subset of rows has the target SNP's value
    and mask zeroed out, and the reconstruction loss is computed only
    on those hidden positions.
    """
    X_tensor = torch.from_numpy(X_norm_filled).to(device)
    mask_tensor = torch.from_numpy(missing_mask.astype(np.float32)).to(device)

    dataset = TensorDataset(X_tensor, mask_tensor)
    loader = DataLoader(dataset, batch_size=cfg.batch_size, shuffle=True)

    n_features = X_norm_filled.shape[1]
    model = AE(n_features, bottleneck=cfg.bottleneck, hidden_dim=cfg.hidden_dim).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=cfg.learning_rate)
    criterion = nn.MSELoss(reduction="none")

    model.train()
    for epoch in range(cfg.n_epochs):
        epoch_loss = 0.0
        for xb, mb in loader:
            optimizer.zero_grad()

            xb_in = xb.clone()
            mb_in = mb.clone()

            hide_rows = torch.rand(xb_in.size(0), device=xb_in.device) < cfg.p_hide

            mb_in[hide_rows, target_j] = 0.0
            xb_in[hide_rows, target_j] = 0.0

            recon = model(xb_in)

            impute_mask = torch.zeros_like(mb_in)
            impute_mask[hide_rows, target_j] = 1.0

            loss_mat = criterion(recon, xb)  # compare against true (unmasked) xb
            masked_loss = (loss_mat * impute_mask).sum() / impute_mask.sum().clamp(min=1.0)

            masked_loss.backward()
            optimizer.step()

            epoch_loss += masked_loss.item()

        epoch_loss /= len(loader)
        if verbose:
            print(f"Epoch {epoch + 1}/{cfg.n_epochs} - impute-loss: {epoch_loss:.4f}")

    return model


def impute_target_snp(
    model: AE,
    X_norm_filled: np.ndarray,
    target_j: int,
    snp_mean: np.ndarray,
    snp_std: np.ndarray,
    device: torch.device,
) -> np.ndarray:
    """
    Run the trained model on all samples with the target SNP zeroed
    out, then undo normalization to return dosages on the original
    0-2 scale.
    """
    X_infer = X_norm_filled.copy()
    X_infer[:, target_j] = 0.0

    X_infer_tensor = torch.from_numpy(X_infer).to(device)

    model.eval()
    with torch.no_grad():
        recon_all = model(X_infer_tensor).cpu().numpy()

    return recon_all * snp_std + snp_mean
