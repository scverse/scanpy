"""Shared docstrings for experimental function parameters."""

from __future__ import annotations

doc_adata = """\
adata
    The annotated data matrix of shape `n_obs` × `n_vars`.
    Rows correspond to cells and columns to genes.
"""

doc_dist_params = """\
theta
    The negative binomial overdispersion parameter `theta` for Pearson residuals.
    Higher values correspond to less overdispersion \
    (`var = mean + mean^2/theta`), and `theta=np.inf` corresponds to a Poisson model.
clip
    Determines if and how residuals are clipped:

    * If `None`, residuals are clipped to the interval \
    `[-sqrt(n_obs), sqrt(n_obs)]`, where `n_obs` is the number of cells in the dataset (default behavior).
    * If any scalar `c`, residuals are clipped to the interval `[-c, c]`. Set \
    `clip=np.inf` for no clipping.
"""

doc_check_values = """\
check_values
    If `True`, checks if counts in selected layer are integers as expected by this
    function, and return a warning if non-integers are found. Otherwise, proceed
    without checking. Setting this to `False` can speed up code for large datasets.
"""

doc_layer = """\
layer
    Layer to use as input instead of `X`. If `None`, `X` is used.
"""

doc_subset = """\
subset
    Inplace subset to highly-variable genes if `True` otherwise merely indicate
    highly variable genes.
"""

doc_genes_batch_chunk = """\
n_top_genes
    Number of highly-variable genes to keep. Mandatory if `flavor='seurat_v3'` or
    `flavor='pearson_residuals'`.
batch_key
    If specified, highly-variable genes are selected within each batch separately
    and merged. This simple process avoids the selection of batch-specific genes
    and acts as a lightweight batch correction method. Genes are first sorted by
    how many batches they are a HVG. If `flavor='pearson_residuals'`, ties are
    broken by the median rank (across batches) based on within-batch residual
    variance.
chunksize
    If `flavor='pearson_residuals'`, this dertermines how many genes are processed at
    once while computing the residual variance. Choosing a smaller value will reduce
    the required memory.
"""

doc_pca_chunk = """\
n_comps
    Number of principal components to compute in the PCA step.
random_state
    Random seed for setting the initial states for the optimization in the PCA step.
kwargs_pca
    Dictionary of further keyword arguments passed on to `scanpy.pp.pca()`.
"""

doc_inplace = """\
inplace
    If `True`, update `adata` with results. Otherwise, return results. See below for
    details of what is returned.
"""

doc_copy = """\
copy
    If `True`, the function runs on a copy of the input object and returns the
    modified copy. Otherwise, the input object is modified direcly. Not compatible
    with `inplace=False`.
"""
