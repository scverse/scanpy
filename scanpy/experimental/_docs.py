"""Shared docstrings for experimental function parameters.
"""

doc_adata = """\
adata
    The annotated data matrix of shape `n_obs` × `n_vars`.
    Rows correspond to cells and columns to genes.
"""

doc_dist_params = """\
theta
    The negative binomial overdispersion parameter theta for Pearson residuals.
    Higher values correspond to less overdispersion (var = mean + mean^2/theta),
    and `theta=np.Inf` corresponds to a Poisson model.
clip
    Determines if and how residuals are clipped:

    * If `None`, residuals are clipped to the interval [-sqrt(n), sqrt(n)], \
    where n is the number of cells in the dataset (default behavior).
    * If any scalar c, residuals are clipped to the interval [-c, c]. Set \
    `clip=np.Inf` for no clipping.
"""

doc_layer = """\
check_values
    Check if counts in selected layer are integers. A Warning is returned if set to
    True.
layer
    Layer to normalize instead of `X`. If `None`, `X` is normalized.
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

doc_inplace = """\
inplace
    Whether to update `adata` or return dictionary with normalized copies
    of `adata.X` and `adata.layers`.
"""

doc_copy = """\
copy
    Whether to modify copied input object. Not compatible with `inplace=False`.
"""
