from typing import Optional, Union

import numpy as np
from scipy.sparse import issparse, spmatrix
from scipy.sparse.linalg import LinearOperator, svds
from sklearn.utils import check_array, check_random_state
from sklearn.utils.extmath import svd_flip

from anndata import AnnData

from .. import logging as logg
from ._utils import _get_mean_var
from .._utils import AnyRandom
from .. import settings


def pca(
    data: Union[AnnData, np.ndarray, spmatrix],
    n_comps: Optional[int] = None,
    zero_center: Optional[bool] = True,
    svd_solver: str = 'arpack',
    random_state: AnyRandom = 0,
    return_info: bool = False,
    use_highly_variable: Optional[bool] = None,
    dtype: str = 'float32',
    copy: bool = False,
    chunked: bool = False,
    chunk_size: Optional[int] = None,
) -> Union[AnnData, np.ndarray, spmatrix]:
    """\
    Principal component analysis [Pedregosa11]_.

    Computes PCA coordinates, loadings and variance decomposition.
    Uses the implementation of *scikit-learn* [Pedregosa11]_.

    .. versionchanged:: 1.5.0

        In previous versions, computing a PCA on a sparse matrix would make a dense copy of
        the array for mean centering.
        As of scanpy 1.5.0, mean centering is implicit.
        While results are extremely similar, they are not exactly the same.
        If you would like to reproduce the old results, pass a dense array.

    Parameters
    ----------
    data
        The (annotated) data matrix of shape `n_obs` × `n_vars`.
        Rows correspond to cells and columns to genes.
    n_comps
        Number of principal components to compute. Defaults to 50, or 1 - minimum
        dimension size of selected representation.
    zero_center
        If `True`, compute standard PCA from covariance matrix.
        If `False`, omit zero-centering variables
        (uses :class:`~sklearn.decomposition.TruncatedSVD`),
        which allows to handle sparse input efficiently.
        Passing `None` decides automatically based on sparseness of the data.
    svd_solver
        SVD solver to use:

        `'arpack'` (the default)
          for the ARPACK wrapper in SciPy (:func:`~scipy.sparse.linalg.svds`)
        `'randomized'`
          for the randomized algorithm due to Halko (2009).
        `'auto'`
          chooses automatically depending on the size of the problem.
        `'lobpcg'`
          An alternative SciPy solver.

        .. versionchanged:: 1.4.5
           Default value changed from `'auto'` to `'arpack'`.

        Efficient computation of the principal components of a sparse matrix
        currently only works with the `'arpack`' or `'lobpcg'` solvers.

    random_state
        Change to use different initial states for the optimization.
    return_info
        Only relevant when not passing an :class:`~anndata.AnnData`:
        see “**Returns**”.
    use_highly_variable
        Whether to use highly variable genes only, stored in
        `.var['highly_variable']`.
        By default uses them if they have been determined beforehand.
    dtype
        Numpy data type string to which to convert the result.
    copy
        If an :class:`~anndata.AnnData` is passed, determines whether a copy
        is returned. Is ignored otherwise.
    chunked
        If `True`, perform an incremental PCA on segments of `chunk_size`.
        The incremental PCA automatically zero centers and ignores settings of
        `random_seed` and `svd_solver`. If `False`, perform a full PCA.
    chunk_size
        Number of observations to include in each chunk.
        Required if `chunked=True` was passed.

    Returns
    -------
    X_pca : :class:`~scipy.sparse.spmatrix`, :class:`~numpy.ndarray`
        If `data` is array-like and `return_info=False` was passed,
        this function only returns `X_pca`…
    adata : anndata.AnnData
        …otherwise if `copy=True` it returns or else adds fields to `adata`:

        `.obsm['X_pca']`
             PCA representation of data.
        `.varm['PCs']`
             The principal components containing the loadings.
        `.uns['pca']['variance_ratio']`
             Ratio of explained variance.
        `.uns['pca']['variance']`
             Explained variance, equivalent to the eigenvalues of the
             covariance matrix.
    """
    logg_start = logg.info('computing PCA')

    # chunked calculation is not randomized, anyways
    if svd_solver in {'auto', 'randomized'} and not chunked:
        logg.info(
            'Note that scikit-learn\'s randomized PCA might not be exactly '
            'reproducible across different computational platforms. For exact '
            'reproducibility, choose `svd_solver=\'arpack\'.`'
        )
    data_is_AnnData = isinstance(data, AnnData)
    if data_is_AnnData:
        adata = data.copy() if copy else data
    else:
        adata = AnnData(data, dtype=data.dtype)

    if use_highly_variable is True and 'highly_variable' not in adata.var.keys():
        raise ValueError(
            'Did not find adata.var[\'highly_variable\']. '
            'Either your data already only consists of highly-variable genes '
            'or consider running `pp.highly_variable_genes` first.'
        )
    if use_highly_variable is None:
        use_highly_variable = True if 'highly_variable' in adata.var.keys() else False
    if use_highly_variable:
        logg.info('    on highly variable genes')
    adata_comp = (
        adata[:, adata.var['highly_variable']] if use_highly_variable else adata
    )

    if n_comps is None:
        min_dim = min(adata_comp.n_vars, adata_comp.n_obs)
        if settings.N_PCS >= min_dim:
            n_comps = min_dim - 1
        else:
            n_comps = settings.N_PCS

    logg.info(f'    with n_comps={n_comps}')

    random_state = check_random_state(random_state)

    X = adata_comp.X

    if chunked:
        if not zero_center or random_state or svd_solver != 'arpack':
            logg.debug('Ignoring zero_center, random_state, svd_solver')

        from sklearn.decomposition import IncrementalPCA

        X_pca = np.zeros((X.shape[0], n_comps), X.dtype)

        pca_ = IncrementalPCA(n_components=n_comps)

        for chunk, _, _ in adata_comp.chunked_X(chunk_size):
            chunk = chunk.toarray() if issparse(chunk) else chunk
            pca_.partial_fit(chunk)

        for chunk, start, end in adata_comp.chunked_X(chunk_size):
            chunk = chunk.toarray() if issparse(chunk) else chunk
            X_pca[start:end] = pca_.transform(chunk)
    elif (not issparse(X) or svd_solver == "randomized") and zero_center:
        from sklearn.decomposition import PCA

        if issparse(X) and svd_solver == "randomized":
            # This  is for backwards compat. Better behaviour would be to either error or use arpack.
            logg.warning(
                "svd_solver 'randomized' does not work with sparse input. Densifying the array. "
                "This may take a very large amount of memory."
            )
            X = X.toarray()
        pca_ = PCA(
            n_components=n_comps, svd_solver=svd_solver, random_state=random_state
        )
        X_pca = pca_.fit_transform(X)
    elif issparse(X) and zero_center:
        from sklearn.decomposition import PCA

        if svd_solver == "auto":
            svd_solver = "arpack"
        if svd_solver not in {'lobpcg', 'arpack'}:
            raise ValueError(
                'svd_solver: {svd_solver} can not be used with sparse input.\n'
                'Use "arpack" (the default) or "lobpcg" instead.'
            )

        output = _pca_with_sparse(
            X, n_comps, solver=svd_solver, random_state=random_state
        )
        # this is just a wrapper for the results
        X_pca = output['X_pca']
        pca_ = PCA(n_components=n_comps, svd_solver=svd_solver)
        pca_.components_ = output['components']
        pca_.explained_variance_ = output['variance']
        pca_.explained_variance_ratio_ = output['variance_ratio']
    elif not zero_center:
        from sklearn.decomposition import TruncatedSVD

        logg.debug(
            '    without zero-centering: \n'
            '    the explained variance does not correspond to the exact statistical defintion\n'
            '    the first component, e.g., might be heavily influenced by different means\n'
            '    the following components often resemble the exact PCA very closely'
        )
        pca_ = TruncatedSVD(
            n_components=n_comps, random_state=random_state, algorithm=svd_solver
        )
        X_pca = pca_.fit_transform(X)
    else:
        raise Exception("This shouldn't happen. Please open a bug report.")

    if X_pca.dtype.descr != np.dtype(dtype).descr:
        X_pca = X_pca.astype(dtype)

    if data_is_AnnData:
        adata.obsm['X_pca'] = X_pca
        adata.uns['pca'] = {}
        adata.uns['pca']['params'] = {
            'zero_center': zero_center,
            'use_highly_variable': use_highly_variable,
        }
        if use_highly_variable:
            adata.varm['PCs'] = np.zeros(shape=(adata.n_vars, n_comps))
            adata.varm['PCs'][adata.var['highly_variable']] = pca_.components_.T
        else:
            adata.varm['PCs'] = pca_.components_.T
        adata.uns['pca']['variance'] = pca_.explained_variance_
        adata.uns['pca']['variance_ratio'] = pca_.explained_variance_ratio_
        logg.info('    finished', time=logg_start)
        logg.debug(
            'and added\n'
            '    \'X_pca\', the PCA coordinates (adata.obs)\n'
            '    \'PC1\', \'PC2\', ..., the loadings (adata.var)\n'
            '    \'pca_variance\', the variance / eigenvalues (adata.uns)\n'
            '    \'pca_variance_ratio\', the variance ratio (adata.uns)'
        )
        return adata if copy else None
    else:
        logg.info('    finished', time=logg_start)
        if return_info:
            return (
                X_pca,
                pca_.components_,
                pca_.explained_variance_ratio_,
                pca_.explained_variance_,
            )
        else:
            return X_pca


def _pca_with_sparse(X, npcs, solver='arpack', mu=None, random_state=None):
    random_state = check_random_state(random_state)
    np.random.set_state(random_state.get_state())
    random_init = np.random.rand(np.min(X.shape))
    X = check_array(X, accept_sparse=['csr', 'csc'])

    if mu is None:
        mu = X.mean(0).A.flatten()[None, :]
    mdot = mu.dot
    mmat = mdot
    mhdot = mu.T.dot
    mhmat = mu.T.dot
    Xdot = X.dot
    Xmat = Xdot
    XHdot = X.T.conj().dot
    XHmat = XHdot
    ones = np.ones(X.shape[0])[None, :].dot

    def matvec(x):
        return Xdot(x) - mdot(x)

    def matmat(x):
        return Xmat(x) - mmat(x)

    def rmatvec(x):
        return XHdot(x) - mhdot(ones(x))

    def rmatmat(x):
        return XHmat(x) - mhmat(ones(x))

    XL = LinearOperator(
        matvec=matvec,
        dtype=X.dtype,
        matmat=matmat,
        shape=X.shape,
        rmatvec=rmatvec,
        rmatmat=rmatmat,
    )

    u, s, v = svds(XL, solver=solver, k=npcs, v0=random_init)
    u, v = svd_flip(u, v)
    idx = np.argsort(-s)
    v = v[idx, :]

    X_pca = (u * s)[:, idx]
    ev = s[idx] ** 2 / (X.shape[0] - 1)

    total_var = _get_mean_var(X)[1].sum()
    ev_ratio = ev / total_var

    output = {
        'X_pca': X_pca,
        'variance': ev,
        'variance_ratio': ev_ratio,
        'components': v,
    }
    return output
