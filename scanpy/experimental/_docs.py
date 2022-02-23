"""Shared docstrings for experimental function parameters.
"""

doc_adata = """\
adata
    The annotated data matrix of shape `n_obs` × `n_vars`.
    Rows correspond to cells and columns to genes.
"""

doc_norm_params = """\
theta
    The negative binomial overdispersion parameter theta for Pearson residuals.
    Higher values correspond to less overdispersion (var = mean + mean^2/theta),
    and `theta=np.Inf` corresponds to a Poisson model.
clip
    Determines if and how residuals are clipped:

    * If `None`, residuals are clipped to the interval [-sqrt(n), sqrt(n)], \
    where n is the number of cells in the dataset (default behavior).
    * If any scalar c, residuals are clipped to the interval [-c, c]. Set \
    `clip=np.Inf` for no clipping.0

"""

doc_layer_copy = """\
check_values
    Check if counts in selected layer are integers. A Warning is returned if set to
    True.
layer
    Layer to normalize instead of `X`. If `None`, `X` is normalized.
copy
    Whether to modify copied input object. Not compatible with `inplace=False`.
inplace
    Whether to update `adata` or return dictionary with normalized copies
    of `adata.X` and `adata.layers`.
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
