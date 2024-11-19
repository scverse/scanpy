from __future__ import annotations

import warnings
from typing import TYPE_CHECKING, Literal, overload
from warnings import warn

import anndata as ad
import numpy as np
from anndata import AnnData
from packaging.version import Version
from scipy.sparse import issparse
from sklearn.utils import check_random_state

from ... import logging as logg
from ..._compat import DaskArray, pkg_version
from ..._settings import settings
from ..._utils import _doc_params, _empty, get_literal_vals, is_backed_type
from ...get import _check_mask, _get_obs_rep
from .._docs import doc_mask_var_hvg
from ._compat import _pca_compat_sparse

if TYPE_CHECKING:
    from collections.abc import Container
    from collections.abc import Set as AbstractSet
    from typing import LiteralString, TypeVar

    import dask_ml.decomposition as dmld
    import sklearn.decomposition as skld
    from numpy.typing import DTypeLike, NDArray
    from scipy import sparse
    from scipy.sparse import spmatrix

    from ..._compat import _LegacyRandom
    from ..._utils import Empty

    CSMatrix = sparse.csr_matrix | sparse.csc_matrix

    MethodDaskML = type[dmld.PCA | dmld.IncrementalPCA | dmld.TruncatedSVD]
    MethodSklearn = type[skld.PCA | skld.TruncatedSVD]

    T = TypeVar("T", bound=LiteralString)
    M = TypeVar("M", bound=LiteralString)


SvdSolvPCADaskML = Literal["auto", "full", "tsqr", "randomized"]
SvdSolvTruncatedSVDDaskML = Literal["tsqr", "randomized"]
SvdSolvDaskML = SvdSolvPCADaskML | SvdSolvTruncatedSVDDaskML

if pkg_version("scikit-learn") >= Version("1.5") or TYPE_CHECKING:
    SvdSolvPCASparseSklearn = Literal["arpack", "covariance_eigh"]
else:
    SvdSolvPCASparseSklearn = Literal["arpack"]
SvdSolvPCADenseSklearn = Literal["auto", "full", "randomized"] | SvdSolvPCASparseSklearn
SvdSolvTruncatedSVDSklearn = Literal["arpack", "randomized"]
SvdSolvSkearn = (
    SvdSolvPCADenseSklearn | SvdSolvPCASparseSklearn | SvdSolvTruncatedSVDSklearn
)

SvdSolvPCACustom = Literal["covariance_eigh"]

SvdSolver = SvdSolvDaskML | SvdSolvSkearn | SvdSolvPCACustom


@_doc_params(
    mask_var_hvg=doc_mask_var_hvg,
)
def pca(
    data: AnnData | np.ndarray | spmatrix,
    n_comps: int | None = None,
    *,
    layer: str | None = None,
    zero_center: bool | None = True,
    svd_solver: SvdSolver | None = None,
    random_state: _LegacyRandom = 0,
    return_info: bool = False,
    mask_var: NDArray[np.bool_] | str | None | Empty = _empty,
    use_highly_variable: bool | None = None,
    dtype: DTypeLike = "float32",
    chunked: bool = False,
    chunk_size: int | None = None,
    key_added: str | None = None,
    copy: bool = False,
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
            For sparse *dask* arrays, will use `'covariance_eigh'`.
            If *scikit-learn* :class:`~sklearn.decomposition.PCA` is used, will give `'arpack'`,
            if *scikit-learn* :class:`~sklearn.decomposition.TruncatedSVD` is used, will give `'randomized'`,
            if *dask-ml* :class:`~dask_ml.decomposition.PCA` or :class:`~dask_ml.decomposition.IncrementalPCA` is used, will give `'auto'`,
            if *dask-ml* :class:`~dask_ml.decomposition.TruncatedSVD` is used, will give `'tsqr'`
        `'arpack'`
            for the ARPACK wrapper in SciPy (:func:`~scipy.sparse.linalg.svds`)
            Not available with *dask* arrays.
        `'covariance_eigh'`
            Classic eigendecomposition of the covariance matrix, suited for tall-and-skinny matrices.
            With dask, array must be CSR and chunked as (N, adata.shape[1]).
        `'randomized'`
            for the randomized algorithm due to Halko (2009). For *dask* arrays,
            this will use :func:`~dask.array.linalg.svd_compressed`.
        `'auto'`
            chooses automatically depending on the size of the problem.
        `'tsqr'`
            Only available with dense *dask* arrays. "tsqr"
            algorithm from Benson et. al. (2013).

        .. versionchanged:: 1.9.3
           Default value changed from `'arpack'` to None.
        .. versionchanged:: 1.4.5
           Default value changed from `'auto'` to `'arpack'`.

        Efficient computation of the principal components of a sparse matrix
        currently only works with the `'arpack`' or `'covariance_eigh`' solver.

        If X is a sparse *dask* array, a custom `'covariance_eigh'` solver will be used.
        If X is a dense *dask* array, *dask-ml* classes :class:`~dask_ml.decomposition.PCA`,
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
    key_added
        If not specified, the embedding is stored as
        :attr:`~anndata.AnnData.obsm`\\ `['X_pca']`, the loadings as
        :attr:`~anndata.AnnData.varm`\\ `['PCs']`, and the the parameters in
        :attr:`~anndata.AnnData.uns`\\ `['pca']`.
        If specified, the embedding is stored as
        :attr:`~anndata.AnnData.obsm`\\ ``[key_added]``, the loadings as
        :attr:`~anndata.AnnData.varm`\\ ``[key_added]``, and the the parameters in
        :attr:`~anndata.AnnData.uns`\\ ``[key_added]``.
    copy
        If an :class:`~anndata.AnnData` is passed, determines whether a copy
        is returned. Is ignored otherwise.

    Returns
    -------
    If `data` is array-like and `return_info=False` was passed,
    this function returns the PCA representation of `data` as an
    array of the same type as the input array.

    Otherwise, it returns `None` if `copy=False`, else an updated `AnnData` object.
    Sets the following fields:

    `.obsm['X_pca' | key_added]` : :class:`~scipy.sparse.spmatrix` | :class:`~numpy.ndarray` (shape `(adata.n_obs, n_comps)`)
        PCA representation of data.
    `.varm['PCs' | key_added]` : :class:`~numpy.ndarray` (shape `(adata.n_vars, n_comps)`)
        The principal components containing the loadings.
    `.uns['pca' | key_added]['variance_ratio']` : :class:`~numpy.ndarray` (shape `(n_comps,)`)
        Ratio of explained variance.
    `.uns['pca' | key_added]['variance']` : :class:`~numpy.ndarray` (shape `(n_comps,)`)
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
            "reproducibility, choose `svd_solver='arpack'`."
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
        n_comps = min_dim - 1 if min_dim <= settings.N_PCS else settings.N_PCS

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

    # check_random_state returns a numpy RandomState when passed an int but
    # dask needs an int for random state
    if not isinstance(X, DaskArray):
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
        if isinstance(X, DaskArray):
            from dask.array import zeros
            from dask_ml.decomposition import IncrementalPCA

            incremental_pca_kwargs["svd_solver"] = _handle_dask_ml_args(
                svd_solver, IncrementalPCA
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
    elif zero_center:
        if issparse(X) and (
            pkg_version("scikit-learn") < Version("1.4") or svd_solver == "lobpcg"
        ):
            if svd_solver not in (
                {"lobpcg"} | get_literal_vals(SvdSolvPCASparseSklearn)
            ):
                if svd_solver is not None:
                    msg = (
                        f"Ignoring {svd_solver=} and using 'arpack', "
                        "sparse PCA with sklearn < 1.4 only supports 'lobpcg' and 'arpack'."
                    )
                    warnings.warn(msg)
                svd_solver = "arpack"
            elif svd_solver == "lobpcg":
                msg = (
                    f"{svd_solver=} for sparse relies on legacy code and will not be supported in the future. "
                    "Also the lobpcg solver has been observed to be inaccurate. Please use 'arpack' instead."
                )
                warnings.warn(msg, FutureWarning)
            X_pca, pca_ = _pca_compat_sparse(
                X, n_comps, solver=svd_solver, random_state=random_state
            )
        else:
            if not isinstance(X, DaskArray):
                from sklearn.decomposition import PCA

                svd_solver = _handle_sklearn_args(svd_solver, PCA, sparse=issparse(X))
                pca_ = PCA(
                    n_components=n_comps,
                    svd_solver=svd_solver,
                    random_state=random_state,
                )
            elif issparse(X._meta):
                from ._dask_sparse import PCASparseDask

                if random_state != 0:
                    msg = f"Ignoring {random_state=} when using a sparse dask array"
                    warnings.warn(msg)
                if svd_solver not in {None, "covariance_eigh"}:
                    msg = f"Ignoring {svd_solver=} when using a sparse dask array"
                    warnings.warn(msg)
                pca_ = PCASparseDask(n_components=n_comps)
            else:
                from dask_ml.decomposition import PCA

                svd_solver = _handle_dask_ml_args(svd_solver, PCA)
                pca_ = PCA(
                    n_components=n_comps,
                    svd_solver=svd_solver,
                    random_state=random_state,
                )
            X_pca = pca_.fit_transform(X)
    else:
        if isinstance(X, DaskArray):
            if issparse(X._meta):
                msg = "Dask sparse arrays do not support zero-centering (yet)"
                raise TypeError(msg)
            from dask_ml.decomposition import TruncatedSVD

            svd_solver = _handle_dask_ml_args(svd_solver, TruncatedSVD)
        else:
            from sklearn.decomposition import TruncatedSVD

            svd_solver = _handle_sklearn_args(svd_solver, TruncatedSVD)

        logg.debug(
            "    without zero-centering: \n"
            "    the explained variance does not correspond to the exact statistical definition\n"
            "    the first component, e.g., might be heavily influenced by different means\n"
            "    the following components often resemble the exact PCA very closely"
        )
        pca_ = TruncatedSVD(
            n_components=n_comps, random_state=random_state, algorithm=svd_solver
        )
        X_pca = pca_.fit_transform(X)

    if X_pca.dtype.descr != np.dtype(dtype).descr:
        X_pca = X_pca.astype(dtype)

    if data_is_AnnData:
        key_obsm, key_varm, key_uns = (
            ("X_pca", "PCs", "pca") if key_added is None else [key_added] * 3
        )
        adata.obsm[key_obsm] = X_pca

        if mask_var is not None:
            adata.varm[key_varm] = np.zeros(shape=(adata.n_vars, n_comps))
            adata.varm[key_varm][mask_var] = pca_.components_.T
        else:
            adata.varm[key_varm] = pca_.components_.T

        params = dict(
            zero_center=zero_center,
            use_highly_variable=mask_var_param == "highly_variable",
            mask_var=mask_var_param,
        )
        if layer is not None:
            params["layer"] = layer
        adata.uns[key_uns] = dict(
            params=params,
            variance=pca_.explained_variance_,
            variance_ratio=pca_.explained_variance_ratio_,
        )

        logg.info("    finished", time=logg_start)
        logg.debug(
            "and added\n"
            f"    {key_obsm!r}, the PCA coordinates (adata.obs)\n"
            f"    {key_varm!r}, the loadings (adata.varm)\n"
            f"    'pca_variance', the variance / eigenvalues (adata.uns[{key_uns!r}])\n"
            f"    'pca_variance_ratio', the variance ratio (adata.uns[{key_uns!r}])"
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
        and "highly_variable" in adata.var.columns
    ):
        mask_var = "highly_variable"

    # Without highly variable genes, we don’t use a mask by default
    if mask_var is _empty or mask_var is None:
        return None, None
    return mask_var, _check_mask(adata, mask_var, "var")


@overload
def _handle_dask_ml_args(
    svd_solver: str | None, method: type[dmld.PCA | dmld.IncrementalPCA]
) -> SvdSolvPCADaskML: ...
@overload
def _handle_dask_ml_args(
    svd_solver: str | None, method: type[dmld.TruncatedSVD]
) -> SvdSolvTruncatedSVDDaskML: ...
def _handle_dask_ml_args(svd_solver: str | None, method: MethodDaskML) -> SvdSolvDaskML:
    import dask_ml.decomposition as dmld

    args: AbstractSet[SvdSolvDaskML]
    default: SvdSolvDaskML
    match method:
        case dmld.PCA | dmld.IncrementalPCA:
            args = get_literal_vals(SvdSolvPCADaskML)
            default = "auto"
        case dmld.TruncatedSVD:
            args = get_literal_vals(SvdSolvTruncatedSVDDaskML)
            default = "tsqr"
        case _:
            msg = f"Unknown {method=} in _handle_dask_ml_args"
            raise ValueError(msg)
    return _handle_x_args(svd_solver, method, args, default)


@overload
def _handle_sklearn_args(
    svd_solver: str | None, method: type[skld.TruncatedSVD], *, sparse: None = None
) -> SvdSolvTruncatedSVDSklearn: ...
@overload
def _handle_sklearn_args(
    svd_solver: str | None, method: type[skld.PCA], *, sparse: Literal[False]
) -> SvdSolvPCADenseSklearn: ...
@overload
def _handle_sklearn_args(
    svd_solver: str | None, method: type[skld.PCA], *, sparse: Literal[True]
) -> SvdSolvPCASparseSklearn: ...
def _handle_sklearn_args(
    svd_solver: str | None, method: MethodSklearn, *, sparse: bool | None = None
) -> SvdSolvSkearn:
    import sklearn.decomposition as skld

    args: AbstractSet[SvdSolvSkearn]
    default: SvdSolvSkearn
    suffix = ""
    match (method, sparse):
        case (skld.TruncatedSVD, None):
            args = get_literal_vals(SvdSolvTruncatedSVDSklearn)
            default = "randomized"
        case (skld.PCA, False):
            args = get_literal_vals(SvdSolvPCADenseSklearn)
            default = "arpack"
        case (skld.PCA, True):
            args = get_literal_vals(SvdSolvPCASparseSklearn)
            default = "arpack"
            suffix = " (with sparse input)"
        case _:
            msg = f"Unknown {method=} ({sparse=}) in _handle_sklearn_args"
            raise ValueError(msg)

    return _handle_x_args(svd_solver, method, args, default, suffix=suffix)


def _handle_x_args(
    svd_solver: str | None,
    method: type,
    args: Container[T],
    default: T,
    *,
    suffix: str = "",
) -> T:
    if svd_solver in args:
        return svd_solver
    if svd_solver is not None:
        msg = (
            f"Ignoring {svd_solver=} and using {default}, "
            f"{method.__module__}.{method.__qualname__}{suffix} only supports {args}."
        )
        warnings.warn(msg)
    return default
