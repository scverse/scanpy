doc_norm_bulk = """\
Normalize total counts per cell.

Normalize each cell by total counts over genes, so that every cell has
the same total count after normalization.

Similar functions are used, for example, by Seurat [Satija15]_, Cell Ranger
[Zheng17]_ or SPRING [Weinreb17]_.

Parameters
----------
adata : :class:`~anndata.AnnData`
    The annotated data matrix of shape `n_obs` × `n_vars`. Rows correspond
    to cells and columns to genes.
target_sum : `float` or `None`, optional (default: `None`)
    If `None`, after normalization, each observation (cell) has a total count
    equal to the median of total counts for observations (cells)
    before normalization.
key_added : `str`, optional (default: `None`)
    Name of the field in `adata.obs` where the total counts per cell are
    stored.
layers : `str` or list of `str`, optional (default: `[]`)
    List of layers to normalize. Set to `'all'` to normalize all layers.
layer_norm : `str` or `None`, optional (default: `None`)
    Specifies how to normalize layers:

    * If `None`, after normalization, for each layer in *layers* each cell
    has a total count equal to the median of the *counts_per_cell* before
    normalization of the layer.
    * If `'after'`, for each layer in *layers* each cell has
    a total count equal to target_sum.
    * If `'X'`, for each layer in *layers* each cell has a total count equal
    to the median of total counts for observations (cells) of `adata.X`
    before normalization.
inplace : `bool`, optional (default: `True`)
    Whether to update `adata` or return dictionary with normalized copies of
    `adata.X` and `adata.layers`.\
"""

doc_norm_quant = """\
quantile : `float`, optional (default: 1)
    Only use genes are less than fraction (specified by *quantile*)
    of the total reads in every cell. So only these genes will sum up
    to *target_sum*.\
"""

doc_norm_return = """\
Returns
-------
Returns dictionary with normalized copies of `adata.X` and `adata.layers`
or updates `adata` with normalized version of the original
`adata.X` and `adata.layers`, depending on `inplace`.\
"""

doc_ex_quant = """\
Examples
--------
>>> adata = AnnData(np.array([[1, 0, 1], [3, 0, 1], [5, 6, 1]]))
>>> sc.pp.normalize_quantile(adata, quantile=0.7)
>>> print(adata.X)
[[1.         0.         1.        ]
 [3.         0.         1.        ]
 [0.71428573 0.85714287 0.14285715]]\
"""

doc_ex_total = """\
Examples
--------
>>> adata = AnnData(np.array([[1, 0], [3, 0], [5, 6]]))
>>> print(adata.X.sum(axis=1))
[  1.   3.  11.]
>>> sc.pp.normalize_total(adata, key_added='n_counts')
>>> print(adata.obs)
>>> print(adata.X.sum(axis=1))
   n_counts
0       1.0
1       3.0
2      11.0
[ 3.  3.  3.]
>>> sc.pp.normalize_total(adata, target_sum=1,
>>>                       key_added='n_counts2')
>>> print(adata.obs)
>>> print(adata.X.sum(axis=1))
   n_counts  n_counts2
0       1.0        3.0
1       3.0        3.0
2      11.0        3.0
[ 1.  1.  1.]\
"""
