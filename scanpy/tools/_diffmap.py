from ._dpt import _diffmap


def diffmap(adata, n_comps=15, copy=False):
    """Diffusion Maps [Coifman05]_ [Haghverdi15]_ [Wolf18]_.

    Diffusion maps [Coifman05]_ has been proposed for visualizing single-cell
    data by [Haghverdi15]_. The tool uses the adapted Gaussian kernel suggested
    by [Haghverdi16]_ in the implementation of [Wolf18]_.

    The width ("sigma") of the connectivity kernel is implicitly determined by
    the number of neighbors used to compute the single-cell graph in
    :func:`~scanpy.api.neighbors`. To reproduce the original implementation
    using a Gaussian kernel, use `method=='gauss'` in
    :func:`~scanpy.api.neighbors`. To use an exponential kernel, use the default
    `method=='umap'`. Differences between these options shouldn't usually be
    dramatic.

    Parameters
    ----------
    adata : :class:`~anndata.AnnData`
        Annotated data matrix.
    n_comps : `int`, optional (default: 15)
        The number of dimensions of the representation.
    copy : `bool` (default: `False`)
        Return a copy instead of writing to adata.

    Returns
    -------
    Depending on `copy`, returns or updates `adata` with the following fields.

    **X_diffmap** : :class:`numpy.ndarray` (`adata.obsm`)
        Diffusion map representation of data, which is the right eigen basis of
        the transition matrix with eigenvectors as columns.
    **diffmap_evals** : :class:`numpy.ndarray` (`adata.uns`)
        Array of size (number of eigen vectors). Eigenvalues of transition matrix.
    """
    if 'neighbors' not in adata.uns:
        raise ValueError(
            'You need to run `pp.neighbors` first to compute a neighborhood graph.')
    if n_comps <= 2:
        raise ValueError(
            'Provide any value greater than 2 for `n_comps`. ')
    adata = adata.copy() if copy else adata
    _diffmap(adata, n_comps=n_comps)
    return adata if copy else None
