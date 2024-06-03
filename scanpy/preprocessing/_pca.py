from __future__ import annotations

import warnings
from typing import TYPE_CHECKING
from warnings import warn

import anndata as ad
import numpy as np
from anndata import AnnData
from packaging.version import Version
from scipy.sparse import issparse
from scipy.sparse.linalg import LinearOperator, svds
from sklearn.utils import check_array, check_random_state
from sklearn.utils.extmath import svd_flip

from .. import logging as logg
from .._compat import DaskArray, pkg_version
from .._settings import settings
from .._utils import _doc_params, _empty, is_backed_type
from ..get import _check_mask, _get_obs_rep
from ._docs import doc_mask_var_hvg
from ._utils import _get_mean_var

if TYPE_CHECKING:
    from numpy.typing import DTypeLike, NDArray
    from scipy.sparse import spmatrix

    from .._utils import AnyRandom, Empty


@_doc_params(
    mask_var_hvg=doc_mask_var_hvg,
)
def pca(
    data: AnnData | np.ndarray | spmatrix,
    n_comps: int | None = None,
    *,
    layer: str | None = None,
    zero_center: bool | None = True,
    svd_solver: str | None = None,
    random_state: AnyRandom = 0,
    return_info: bool = False,
    mask_var: NDArray[np.bool_] | str | None | Empty = _empty,
    use_highly_variable: bool | None = None,
    dtype: DTypeLike = "float32",
    copy: bool = False,
    chunked: bool = False,
    chunk_size: int | None = None,
) -> AnnData | np.ndarray | spmatrix | None:
    """\
    Principal component analysis :cite:p:`Pedregosa2011`.

    Computes PCA coordinates, loadings and variance decomposition.
    Uses the implementation of *scikit-learn* :cite:p:`Pedregosa2011`.

    .. versionchanged:: 1.5.0

        In previous versions, computing a PCA on a sparse matrix would make
        a dense copy of the array for mean centering.
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
    layer
        If provided, which element of layers to use for PCA.
    zero_center
        If `True`, compute standard PCA from covariance matrix.
        If `False`, omit zero-centering variables
        (uses *scikit-learn* :class:`~sklearn.decomposition.TruncatedSVD` or
        *dask-ml* :class:`~dask_ml.decomposition.TruncatedSVD`),
        which allows to handle sparse input efficiently.
        Passing `None` decides automatically based on sparseness of the data.
    svd_solver
        SVD solver to use:

        `None`
            See `chunked` and `zero_center` descriptions to determine which class will be used.
            Depending on the class and the type of X different values for default will be set.
            If *scikit-learn* :class:`~sklearn.decomposition.PCA` is used, will give `'arpack'`,
            if *scikit-learn* :class:`~sklearn.decomposition.TruncatedSVD` is used, will give `'randomized'`,
            if *dask-ml* :class:`~dask_ml.decomposition.PCA` or :class:`~dask_ml.decomposition.IncrementalPCA` is used, will give `'auto'`,
            if *dask-ml* :class:`~dask_ml.decomposition.TruncatedSVD` is used, will give `'tsqr'`
        `'arpack'`
            for the ARPACK wrapper in SciPy (:func:`~scipy.sparse.linalg.svds`)
            Not available with *dask* arrays.
        `'randomized'`
            for the randomized algorithm due to Halko (2009). For *dask* arrays,
            this will use :func:`~dask.array.linalg.svd_compressed`.
        `'auto'`
            chooses automatically depending on the size of the problem.
        `'lobpcg'`
            An alternative SciPy solver. Not available with dask arrays.
        `'tsqr'`
            Only available with *dask* arrays. "tsqr"
            algorithm from Benson et. al. (2013).

        .. versionchanged:: 1.9.3
           Default value changed from `'arpack'` to None.
        .. versionchanged:: 1.4.5
           Default value changed from `'auto'` to `'arpack'`.

        Efficient computation of the principal components of a sparse matrix
        currently only works with the `'arpack`' or `'lobpcg'` solvers.

        If X is a *dask* array, *dask-ml* classes :class:`~dask_ml.decomposition.PCA`,
        :class:`~dask_ml.decomposition.IncrementalPCA`, or
        :class:`~dask_ml.decomposition.TruncatedSVD` will be used.
        Otherwise their *scikit-learn* counterparts :class:`~sklearn.decomposition.PCA`,
        :class:`~sklearn.decomposition.IncrementalPCA`, or
        :class:`~sklearn.decomposition.TruncatedSVD` will be used.
    random_state
        Change to use different initial states for the optimization.
    return_info
        Only relevant when not passing an :class:`~anndata.AnnData`:
        see “Returns”.
    {mask_var_hvg}
    layer
        Layer of `adata` to use as expression values.
    dtype
        Numpy data type string to which to convert the result.
    copy
        If an :class:`~anndata.AnnData` is passed, determines whether a copy
        is returned. Is ignored otherwise.
    chunked
        If `True`, perform an incremental PCA on segments of `chunk_size`.
        The incremental PCA automatically zero centers and ignores settings of
        `random_seed` and `svd_solver`. Uses sklearn :class:`~sklearn.decomposition.IncrementalPCA` or
        *dask-ml* :class:`~dask_ml.decomposition.IncrementalPCA`. If `False`, perform a full PCA and
        use sklearn :class:`~sklearn.decomposition.PCA` or
        *dask-ml* :class:`~dask_ml.decomposition.PCA`
    chunk_size
        Number of observations to include in each chunk.
        Required if `chunked=True` was passed.

    Returns
    -------
    If `data` is array-like and `return_info=False` was passed,
    this function returns the PCA representation of `data` as an
    array of the same type as the input array.

    Otherwise, it returns `None` if `copy=False`, else an updated `AnnData` object.
    Sets the following fields:

    `.obsm['X_pca']` : :class:`~scipy.sparse.spmatrix` | :class:`~numpy.ndarray` (shape `(adata.n_obs, n_comps)`)
        PCA representation of data.
    `.varm['PCs']` : :class:`~numpy.ndarray` (shape `(adata.n_vars, n_comps)`)
        The principal components containing the loadings.
    `.uns['pca']['variance_ratio']` : :class:`~numpy.ndarray` (shape `(n_comps,)`)
        Ratio of explained variance.
    `.uns['pca']['variance']` : :class:`~numpy.ndarray` (shape `(n_comps,)`)
        Explained variance, equivalent to the eigenvalues of the
        covariance matrix.
    """
    logg_start = logg.info("computing PCA")
    if layer is not None and chunked:
        # Current chunking implementation relies on pca being called on X
        raise NotImplementedError("Cannot use `layer` and `chunked` at the same time.")

    # chunked calculation is not randomized, anyways
    if svd_solver in {"auto", "randomized"} and not chunked:
        logg.info(
            "Note that scikit-learn's randomized PCA might not be exactly "
            "reproducible across different computational platforms. For exact "
            "reproducibility, choose `svd_solver='arpack'.`"
        )
    data_is_AnnData = isinstance(data, AnnData)
    if data_is_AnnData:
        if layer is None and not chunked and is_backed_type(data.X):
            raise NotImplementedError(
                f"PCA is not implemented for matrices of type {type(data.X)} with chunked as False"
            )
        adata = data.copy() if copy else data
    else:
        if pkg_version("anndata") < Version("0.8.0rc1"):
            adata = AnnData(data, dtype=data.dtype)
        else:
            adata = AnnData(data)

    # Unify new mask argument and deprecated use_highly_varible argument
    mask_var_param, mask_var = _handle_mask_var(adata, mask_var, use_highly_variable)
    del use_highly_variable
    adata_comp = adata[:, mask_var] if mask_var is not None else adata

    if n_comps is None:
        min_dim = min(adata_comp.n_vars, adata_comp.n_obs)
        if settings.N_PCS >= min_dim:
            n_comps = min_dim - 1
        else:
            n_comps = settings.N_PCS

    logg.info(f"    with n_comps={n_comps}")

    X = _get_obs_rep(adata_comp, layer=layer)
    if is_backed_type(X) and layer is not None:
        raise NotImplementedError(
            f"PCA is not implemented for matrices of type {type(X)} from layers"
        )
    # See: https://github.com/scverse/scanpy/pull/2816#issuecomment-1932650529
    if (
        Version(ad.__version__) < Version("0.9")
        and mask_var is not None
        and isinstance(X, np.ndarray)
    ):
        warnings.warn(
            "When using a mask parameter with anndata<0.9 on a dense array, the PCA"
            "can have slightly different results due the array being column major "
            "instead of row major.",
            UserWarning,
        )

    is_dask = isinstance(X, DaskArray)

    # check_random_state returns a numpy RandomState when passed an int but
    # dask needs an int for random state
    if not is_dask:
        random_state = check_random_state(random_state)
    elif not isinstance(random_state, int):
        msg = f"random_state needs to be an int, not a {type(random_state).__name__} when passing a dask array"
        raise TypeError(msg)

    if chunked:
        if (
            not zero_center
            or random_state
            or (svd_solver is not None and svd_solver != "arpack")
        ):
            logg.debug("Ignoring zero_center, random_state, svd_solver")

        incremental_pca_kwargs = dict()
        if is_dask:
            from dask.array import zeros
            from dask_ml.decomposition import IncrementalPCA

            incremental_pca_kwargs["svd_solver"] = _handle_dask_ml_args(
                svd_solver, "IncrementalPCA"
            )
        else:
            from numpy import zeros
            from sklearn.decomposition import IncrementalPCA

        X_pca = zeros((X.shape[0], n_comps), X.dtype)

        pca_ = IncrementalPCA(n_components=n_comps, **incremental_pca_kwargs)

        for chunk, _, _ in adata_comp.chunked_X(chunk_size):
            chunk = chunk.toarray() if issparse(chunk) else chunk
            pca_.partial_fit(chunk)

        for chunk, start, end in adata_comp.chunked_X(chunk_size):
            chunk = chunk.toarray() if issparse(chunk) else chunk
            X_pca[start:end] = pca_.transform(chunk)
    elif (not issparse(X) or svd_solver == "randomized") and zero_center:
        if is_dask:
            from dask_ml.decomposition import PCA

            svd_solver = _handle_dask_ml_args(svd_solver, "PCA")
        else:
            from sklearn.decomposition import PCA

            svd_solver = _handle_sklearn_args(svd_solver, "PCA")

        if issparse(X) and svd_solver == "randomized":
            # This  is for backwards compat. Better behaviour would be to either error or use arpack.
            warnings.warn(
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

        svd_solver = _handle_sklearn_args(svd_solver, "PCA (with sparse input)")

        output = _pca_with_sparse(
            X, n_comps, solver=svd_solver, random_state=random_state
        )
        # this is just a wrapper for the results
        X_pca = output["X_pca"]
        pca_ = PCA(
            n_components=n_comps, svd_solver=svd_solver, random_state=random_state
        )
        pca_.components_ = output["components"]
        pca_.explained_variance_ = output["variance"]
        pca_.explained_variance_ratio_ = output["variance_ratio"]
    elif not zero_center:
        if is_dask:
            from dask_ml.decomposition import TruncatedSVD

            svd_solver = _handle_dask_ml_args(svd_solver, "TruncatedSVD")
        else:
            from sklearn.decomposition import TruncatedSVD

            svd_solver = _handle_sklearn_args(svd_solver, "TruncatedSVD")

        logg.debug(
            "    without zero-centering: \n"
            "    the explained variance does not correspond to the exact statistical defintion\n"
            "    the first component, e.g., might be heavily influenced by different means\n"
            "    the following components often resemble the exact PCA very closely"
        )
        pca_ = TruncatedSVD(
            n_components=n_comps, random_state=random_state, algorithm=svd_solver
        )
        X_pca = pca_.fit_transform(X)
    else:
        msg = "This shouldn’t happen. Please open a bug report."
        raise AssertionError(msg)

    if X_pca.dtype.descr != np.dtype(dtype).descr:
        X_pca = X_pca.astype(dtype)

    if data_is_AnnData:
        adata.obsm["X_pca"] = X_pca

        if mask_var is not None:
            adata.varm["PCs"] = np.zeros(shape=(adata.n_vars, n_comps))
            adata.varm["PCs"][mask_var] = pca_.components_.T
        else:
            adata.varm["PCs"] = pca_.components_.T

        uns_entry = {
            "params": {
                "zero_center": zero_center,
                "use_highly_variable": mask_var_param == "highly_variable",
                "mask_var": mask_var_param,
            },
            "variance": pca_.explained_variance_,
            "variance_ratio": pca_.explained_variance_ratio_,
        }
        if layer is not None:
            uns_entry["params"]["layer"] = layer
        adata.uns["pca"] = uns_entry

        logg.info("    finished", time=logg_start)
        logg.debug(
            "and added\n"
            "    'X_pca', the PCA coordinates (adata.obs)\n"
            "    'PC1', 'PC2', ..., the loadings (adata.var)\n"
            "    'pca_variance', the variance / eigenvalues (adata.uns)\n"
            "    'pca_variance_ratio', the variance ratio (adata.uns)"
        )
        return adata if copy else None
    else:
        logg.info("    finished", time=logg_start)
        if return_info:
            return (
                X_pca,
                pca_.components_,
                pca_.explained_variance_ratio_,
                pca_.explained_variance_,
            )
        else:
            return X_pca


def _handle_mask_var(
    adata: AnnData,
    mask_var: NDArray[np.bool_] | str | Empty | None,
    use_highly_variable: bool | None,
) -> tuple[np.ndarray | str | None, np.ndarray | None]:
    """\
    Unify new mask argument and deprecated use_highly_varible argument.

    Returns both the normalized mask parameter and the validated mask array.
    """
    # First, verify and possibly warn
    if use_highly_variable is not None:
        hint = (
            'Use_highly_variable=True can be called through mask_var="highly_variable". '
            "Use_highly_variable=False can be called through mask_var=None"
        )
        msg = f"Argument `use_highly_variable` is deprecated, consider using the mask argument. {hint}"
        warn(msg, FutureWarning)
        if mask_var is not _empty:
            msg = f"These arguments are incompatible. {hint}"
            raise ValueError(msg)

    # Handle default case and explicit use_highly_variable=True
    if use_highly_variable or (
        use_highly_variable is None
        and mask_var is _empty
        and "highly_variable" in adata.var.keys()
    ):
        mask_var = "highly_variable"

    # Without highly variable genes, we don’t use a mask by default
    if mask_var is _empty or mask_var is None:
        return None, None
    return mask_var, _check_mask(adata, mask_var, "var")


def _pca_with_sparse(
    X: spmatrix,
    n_pcs: int,
    *,
    solver: str = "arpack",
    mu: NDArray[np.floating] | None = None,
    random_state: AnyRandom = None,
):
    random_state = check_random_state(random_state)
    np.random.set_state(random_state.get_state())
    random_init = np.random.rand(np.min(X.shape))
    X = check_array(X, accept_sparse=["csr", "csc"])

    if mu is None:
        mu = np.asarray(X.mean(0)).flatten()[None, :]
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

    u, s, v = svds(XL, solver=solver, k=n_pcs, v0=random_init)
    # u_based_decision was changed in https://github.com/scikit-learn/scikit-learn/pull/27491
    u, v = svd_flip(
        u, v, u_based_decision=pkg_version("scikit-learn") < Version("1.5.0rc1")
    )
    idx = np.argsort(-s)
    v = v[idx, :]

    X_pca = (u * s)[:, idx]
    ev = s[idx] ** 2 / (X.shape[0] - 1)

    total_var = _get_mean_var(X)[1].sum()
    ev_ratio = ev / total_var

    output = {
        "X_pca": X_pca,
        "variance": ev,
        "variance_ratio": ev_ratio,
        "components": v,
    }
    return output


def _handle_dask_ml_args(svd_solver: str, method: str) -> str:
    method2args = {
        "PCA": {"auto", "full", "tsqr", "randomized"},
        "IncrementalPCA": {"auto", "full", "tsqr", "randomized"},
        "TruncatedSVD": {"tsqr", "randomized"},
    }
    method2default = {
        "PCA": "auto",
        "IncrementalPCA": "auto",
        "TruncatedSVD": "tsqr",
    }

    return _handle_x_args("dask_ml", svd_solver, method, method2args, method2default)


def _handle_sklearn_args(svd_solver: str, method: str) -> str:
    method2args = {
        "PCA": {"auto", "full", "arpack", "randomized"},
        "TruncatedSVD": {"arpack", "randomized"},
        "PCA (with sparse input)": {"lobpcg", "arpack"},
    }
    method2default = {
        "PCA": "arpack",
        "TruncatedSVD": "randomized",
        "PCA (with sparse input)": "arpack",
    }

    return _handle_x_args("sklearn", svd_solver, method, method2args, method2default)


def _handle_x_args(lib, svd_solver, method, method2args, method2default):
    if svd_solver not in method2args[method]:
        if svd_solver is not None:
            warnings.warn(
                f"Ignoring {svd_solver} and using {method2default[method]}, {lib}.decomposition.{method} only supports {method2args[method]}"
            )
        svd_solver = method2default[method]
    return svd_solver
