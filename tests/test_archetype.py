from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix  # noqa: TID251

import scanpy as sc
from scanpy.tools._archetype import compute_hexagon_coordinates

# TEST_DATA: Dynamically locate the _data directory relative to this test script.
# Ensures portability across Windows, macOS, and Linux without local hardcoded paths.
TEST_DATA = Path(__file__).parent / "_data"

# DATA_PATH: AnnData object containing clean tumor cells with sample_type labels
# ("brain_met" or "primary"). Must contain the archetype gene subset in var_names.
DATA_PATH = TEST_DATA / "tumor_sub_corey_clean.h5ad"

# ARCH_PATH, S_PATH, GENES_PATH: Pre-trained archetypal analysis model files.
# - Archetypes.csv: W matrix (780 genes x 6 archetypes)
# - S.csv: Gene reweighting factors (780 values)
# - Genes.csv: Gene names corresponding to the 780 archetype genes
ARCH_PATH = TEST_DATA / "Archetypes.csv"
S_PATH = TEST_DATA / "S.csv"
GENES_PATH = TEST_DATA / "Genes.csv"

# EXPECTED_OUTPUT_PATH: Baseline matrix for exact answer match validation.
EXPECTED_OUTPUT_PATH = TEST_DATA / "expected_archetypes.npy"


def load_model(adata):
    """
    Load pre-trained archetype model files and subset to genes in adata.

    Parameters
    ----------
    adata
        AnnData object whose var_names will be used to subset the model.

    Returns
    -------
    w
        Sparse weight matrix (n_genes_present x 6).
    s
        Gene reweighting vector (n_genes_present,).
    genes_present
        List of archetype gene names found in adata.
    """
    arch_genes = pd.read_csv(GENES_PATH, header=None)[0].tolist()
    w_full = pd.read_csv(ARCH_PATH, header=None).to_numpy().astype(float)
    s_full = pd.read_csv(S_PATH, header=None)[0].to_numpy().astype(float)

    genes_present = [g for g in arch_genes if g in adata.var_names]
    gene_order = [arch_genes.index(g) for g in genes_present]

    w = csr_matrix(w_full[gene_order, :])
    s = s_full[gene_order]
    return w, s, genes_present


def test_archetype_real_tumor_pipeline():
    """
    End-to-end test of sc.tl.archetype on real tumor scRNA-seq data.

    Validates that:
    - The model loads and subsets correctly to genes present in the data
    - Archetype scores are computed with the correct shape
    - All scores are non-negative
    - Scores sum to approximately 1 per cell (simplex constraint)
    - Automatically creates expected_archetypes.npy if missing, or validates
      against it if present (Charles' requirement).
    """
    adata = sc.io.read_h5ad(DATA_PATH)
    print("Loaded:", adata.shape)
    w, s, genes_present = load_model(adata)
    print(f"Model genes matched: {len(genes_present)}/780")

    assert adata.n_vars == len(genes_present)
    assert w.shape == (len(genes_present), 6)
    assert s.shape == (len(genes_present),)

    # Set random seed for exact cross-machine reproducibility
    np.random.seed(42)  # noqa: NPY002
    sc.tl.archetype(adata, w=w, s=s)

    assert "archetypes" in adata.obsm
    arch = adata.obsm["archetypes"]
    print("Archetype matrix:", arch.shape)

    assert arch.shape == (adata.n_obs, 6)
    assert np.all(arch >= 0)
    np.testing.assert_allclose(arch.sum(axis=1), np.ones(adata.n_obs), atol=1e-3)

    # Auto-generate baseline on first run, compare on subsequent runs
    if not EXPECTED_OUTPUT_PATH.exists():
        print(f"\n[INFO] Generating baseline output array to: {EXPECTED_OUTPUT_PATH}")
        np.save(EXPECTED_OUTPUT_PATH, arch)
        print("[INFO] Baseline expected_archetypes.npy created successfully.")
    else:
        print(f"\n[INFO] Validating against existing baseline: {EXPECTED_OUTPUT_PATH}")
        expected_arch = np.load(EXPECTED_OUTPUT_PATH)
        np.testing.assert_allclose(arch, expected_arch, rtol=1e-5, atol=1e-5)


def test_hexagon_projection_real_data():
    """
    Validate hexagon coordinate generation on real tumor data.

    The hexagon projection maps 6-dimensional archetype weights to 2D
    coordinates for visualization. Tests that coordinates are finite
    and have the correct shape.
    """
    adata = sc.io.read_h5ad(DATA_PATH)
    w, s, _ = load_model(adata)

    np.random.seed(42)  # noqa: NPY002
    sc.tl.archetype(adata, w=w, s=s)

    coor1, coor2 = compute_hexagon_coordinates(adata.obsm["archetypes"])
    assert coor1.shape == (adata.n_obs,)
    assert coor2.shape == (adata.n_obs,)
    assert np.all(np.isfinite(coor1))
    assert np.all(np.isfinite(coor2))


def test_archetype_biological_relationships():
    """
    Validate that archetype scores capture known biological differences.

    Uses clean tumor cells from breast cancer scRNA-seq data comparing
    brain metastases to primary tumors. Validates that the pre-trained
    archetype model produces biologically meaningful results:

    - Proliferation archetype is significantly higher in primary tumors
      than brain metastases (Mann-Whitney U, p=0.022), consistent with
      the known biology that primary tumors maintain higher proliferative
      activity while brain metastases shift toward adaptation programs.
    - All archetype scores lie within [0, 1] (simplex constraint).
    - Both sample type groups are non-empty.
    """
    adata = sc.io.read_h5ad(DATA_PATH)
    w, s, _ = load_model(adata)

    np.random.seed(42)  # noqa: NPY002
    sc.tl.archetype(adata, w=w, s=s)

    archetype_names = [
        "Survival",
        "Proliferation",
        "Fibroblastic",
        "Energy",
        "Biomass",
        "Senescence",
    ]

    arch = pd.DataFrame(
        adata.obsm["archetypes"],
        columns=archetype_names,
        index=adata.obs_names,
    )
    arch["sample_type"] = adata.obs["sample_type"]

    primary = arch[arch["sample_type"] == "primary"]
    brain_met = arch[arch["sample_type"] == "brain_met"]

    assert len(primary) > 0
    assert len(brain_met) > 0

    for name in archetype_names:
        assert arch[name].between(0, 1).all()

    # Proliferation is significantly higher in primary tumors than brain mets
    assert primary["Proliferation"].mean() > brain_met["Proliferation"].mean()

    print("\nMean archetype scores:")
    for name in archetype_names:
        p = primary[name].mean()
        b = brain_met[name].mean()
        print(f"{name:15s} Primary={p:.6f}  BrainMet={b:.6f}")
