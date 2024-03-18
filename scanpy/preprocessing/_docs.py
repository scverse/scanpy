"""Shared docstrings for preprocessing function parameters."""

from __future__ import annotations

doc_adata_basic = """\
adata
    Annotated data matrix.\
"""

doc_expr_reps = """\
layer
    If provided, use `adata.layers[layer]` for expression values instead
    of `adata.X`.
use_raw
    If True, use `adata.raw.X` for expression values instead of `adata.X`.\
"""

doc_mask_var_hvg = """\
mask_var
    To run only on a certain set of genes given by a boolean array
    or a string referring to an array in :attr:`~anndata.AnnData.var`.
    By default, uses `.var['highly_variable']` if available, else everything.
use_highly_variable
    Whether to use highly variable genes only, stored in
    `.var['highly_variable']`.
    By default uses them if they have been determined beforehand.

    .. deprecated:: 1.10.0
       Use `mask_var` instead
"""

doc_obs_qc_args = """\
qc_vars
    Keys for boolean columns of `.var` which identify variables you could
    want to control for (e.g. "ERCC" or "mito").
percent_top
    List of ranks (where genes are ranked by expression) at which the cumulative
    proportion of expression will be reported as a percentage. This can be used to
    assess library complexity. Ranks are considered 1-indexed, and if empty or None
    don't calculate.

    E.g. `percent_top=[50]` finds cumulative proportion to the 50th most expressed gene.
"""

doc_qc_metric_naming = """\
expr_type
    Name of kind of values in X.
var_type
    The kind of thing the variables are.\
"""

doc_obs_qc_returns = """\
Observation level metrics include:

`total_{var_type}_by_{expr_type}`
    E.g. "total_genes_by_counts". Number of genes with positive counts in a cell.
`total_{expr_type}`
    E.g. "total_counts". Total number of counts for a cell.
`pct_{expr_type}_in_top_{n}_{var_type}` – for `n` in `percent_top`
    E.g. "pct_counts_in_top_50_genes". Cumulative percentage of counts
    for 50 most expressed genes in a cell.
`total_{expr_type}_{qc_var}` – for `qc_var` in `qc_vars`
    E.g. "total_counts_mito". Total number of counts for variables in
    `qc_vars`.
`pct_{expr_type}_{qc_var}` – for `qc_var` in `qc_vars`
    E.g. "pct_counts_mito". Proportion of total counts for a cell which
    are mitochondrial.\
"""

doc_var_qc_returns = """\
Variable level metrics include:

`total_{expr_type}`
    E.g. "total_counts". Sum of counts for a gene.
`n_genes_by_{expr_type}`
    E.g. "n_genes_by_counts". The number of genes with at least 1 count in a cell. Calculated for all cells.
`mean_{expr_type}`
    E.g. "mean_counts". Mean expression over all cells.
`n_cells_by_{expr_type}`
    E.g. "n_cells_by_counts". Number of cells this expression is
    measured in.
`pct_dropout_by_{expr_type}`
    E.g. "pct_dropout_by_counts". Percentage of cells this feature does
    not appear in.\
"""
