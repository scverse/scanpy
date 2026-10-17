"""Generate PR validation plots for archetype projection."""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import scanpy as sc

# Add repository root to Python path
sys.path.insert(0, str(Path(__file__).parent))

from tests.test_archetype import DATA_PATH, load_model

import scanpy.tools._archetype as arch_mod

# -------------------------------------------------------------------------
# 1. Load Data & Run Archetype Pipeline
# -------------------------------------------------------------------------
print("[INFO] Loading AnnData dataset...")
adata = sc.read_h5ad(DATA_PATH)

print("[INFO] Loading model weights and scaling factors...")
w, s, _ = load_model(adata)

print("[INFO] Computing archetype probabilities...")
np.random.seed(42)  # noqa: NPY002
arch_mod.archetype(adata, w=w, s=s)

# -------------------------------------------------------------------------
# 2. Compute Hexagon Coordinates directly from Archetype Matrix
# -------------------------------------------------------------------------
print("[INFO] Computing hexagon coordinates from archetype matrix...")
arch_matrix = adata.obsm["archetypes"]
res = arch_mod.compute_hexagon_coordinates(arch_matrix)

# Handle tuple return (x_coords, y_coords)
coords = np.column_stack((res[0], res[1])) if isinstance(res, tuple) else res

adata.obsm["X_archetype_hexagon"] = coords

# -------------------------------------------------------------------------
# 3. Extract Archetype Scores & Render Plot
# -------------------------------------------------------------------------
archetype_names = [
    "Survival",
    "Proliferation",
    "Fibroblastic",
    "Energy",
    "Biomass",
    "Senescence",
]
for i, name in enumerate(archetype_names):
    adata.obs[name] = adata.obsm["archetypes"][:, i]

print("[INFO] Rendering hexagon spatial projection plot...")
fig, axes = plt.subplots(2, 3, figsize=(15, 10))
axes = axes.flatten()

for i, name in enumerate(archetype_names):
    ax = axes[i]
    scatter = ax.scatter(
        coords[:, 0],
        coords[:, 1],
        c=adata.obs[name],
        cmap="magma",
        s=3,
        alpha=0.7,
    )
    ax.set_title(f"Archetype: {name}", fontsize=12, fontweight="bold")
    ax.set_aspect("equal")
    ax.axis("off")
    fig.colorbar(scatter, ax=ax, fraction=0.046, pad=0.04)

plt.suptitle("Hexagon Spatial Projection of Tumor Archetypes", fontsize=16)
plt.tight_layout()

output_path = Path("archetype_hexagon_plot.png")
plt.savefig(output_path, dpi=300, bbox_inches="tight")
plt.close()

print(f"[SUCCESS] Visual proof plot saved successfully to: {output_path.resolve()}")
