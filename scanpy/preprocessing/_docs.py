"""Shared docstrings for preprocessing function parameters.
"""

doc_expr_reps = """\
layer
    If provided, use `adata.layers[layer]` for expression values instead
    of `adata.X`.
use_raw
    If True, use `adata.raw.X` for expression values instead of `adata.X`.\
"""

doc_obs_qc_args = """\
qc_vars
    Keys for boolean columns of `.var` which identify variables you could
    want to control for (e.g. "ERCC" or "mito").
percent_top
    Which proportions of top genes to cover. If empty or `None` don't
    calculate. Values are considered 1-indexed, `percent_top=[50]` finds
    cumulative proportion to the 50th most expressed gene.\
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
    E.g. "total_counts_mito". Total number of counts for variabes in
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

doc_adata_basic = """\
adata
    Annotated data matrix.\
"""
