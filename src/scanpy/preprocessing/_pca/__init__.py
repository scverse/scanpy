from __future__ import annotations

from typing import TYPE_CHECKING, Literal, overload

import numpy as np
from anndata import AnnData
from packaging.version import Version

from ... import logging as logg
from ..._compat import CSBase, DaskArray, pkg_version, warn
from ..._settings import settings
from ..._utils import _doc_params, _empty, get_literal_vals, is_backed_type
from ..._utils.random import accepts_legacy_random_state, legacy_random_state
from ...get import _check_mask, _get_obs_rep
from .._docs import doc_mask_var_hvg
from ._compat import _pca_compat_sparse

if TYPE_CHECKING:
    from collections.abc import Container
    from collections.abc import Set as AbstractSet
    from typing import LiteralString

    import dask_ml.decomposition as dmld
    import sklearn.decomposition as skld
    from numpy.typing import DTypeLike, NDArray

    from ..._utils import Empty
    from ..._utils.random import RNGLike, SeedLike


type MethodDaskML = type[dmld.PCA | dmld.IncrementalPCA | dmld.TruncatedSVD]
type MethodSklearn = type[skld.PCA | skld.TruncatedSVD]

type SvdSolvPCADaskML = Literal["auto", "full", "tsqr", "randomized"]
type SvdSolvTruncatedSVDDaskML = Literal["tsqr", "randomized"]
type SvdSolvDaskML = SvdSolvPCADaskML | SvdSolvTruncatedSVDDaskML

if pkg_version("scikit-learn") >= Version("1.5") or TYPE_CHECKING:
    type SvdSolvPCASparseSklearn = Literal["arpack", "covariance_eigh"]
else:
    type SvdSolvPCASparseSklearn = Literal["arpack"]
type SvdSolvPCADenseSklearn = (
    Literal["auto", "full", "randomized"] | SvdSolvPCASparseSklearn
)
type SvdSolvTruncatedSVDSklearn = Literal["arpack", "randomized"]
type SvdSolvSkearn = (
    SvdSolvPCADenseSklearn | SvdSolvPCASparseSklearn | SvdSolvTruncatedSVDSklearn
)

type SvdSolvPCACustom = Literal["covariance_eigh"]
type SvdSolver = SvdSolvDaskML | SvdSolvSkearn | SvdSolvPCACustom


@_doc_params(
    mask_var_hvg=doc_mask_var_hvg,
)
@accepts_legacy_random_state(0)
def pca(  # noqa: PLR0912, PLR0913, PLR0915
    data: AnnData | np.ndarray | CSBase,
    n_comps: int | None = None,
    *,
    layer: str | None = None,
    obsm: str | None = None,
    zero_center: bool = True,
    svd_solver: SvdSolver | None = None,
    chunked: bool = False,
    chunk_size: int | None = None,
    rng: SeedLike | RNGLike | None = None,
    return_info: bool = False,
    mask_var: NDArray[np.bool_] | str | None | Empty = _empty,
    use_highly_variable: bool | None = None,
    dtype: DTypeLike = "float32",
    key_added: str | None = None,
    copy: bool = False,
) -> AnnData | np.ndarray | CSBase | None:
    r"""Principal component analysis :cite:p:`Pedregosa2011`.

    Computes PCA coordinates, loadings and variance decomposition.
    Uses the following implementations (and defaults for `svd_solver`):

    .. list-table::
       :header-rows: 1
       :stub-columns: 1

       - -
         - :class:`~numpy.ndarray`, :class:`~scipy.sparse.spmatrix`, or :class:`~scipy.sparse.sparray`
         - :class:`dask.array.Array`
       - - `chunked=False`, `zero_center=True`
         - sklearn :class:`~sklearn.decomposition.PCA` (`'arpack'`)
         - - *dense*: dask-ml :class:`~dask_ml.decomposition.PCA`\ [#high-mem]_ (`'auto'`)
           - *sparse* or `svd_solver='covariance_eigh'`: custom implementation (`'covariance_eigh'`)
       - - `chunked=False`, `zero_center=False`
         - sklearn :class:`~sklearn.decomposition.TruncatedSVD` (`'randomized'`)
         - dask-ml :class:`~dask_ml.decomposition.TruncatedSVD`\ [#dense-only]_ (`'tsqr'`)
       - - `chunked=True` (`zero_center` ignored)
         - sklearn :class:`~sklearn.decomposition.IncrementalPCA` (`'auto'`)
         - dask-ml :class:`~dask_ml.decomposition.IncrementalPCA`\ [#densifies]_ (`'auto'`)

    .. [#high-mem] Consider `svd_solver='covariance_eigh'` to reduce memory usage (see :issue:`dask/dask-ml#985`).
    .. [#dense-only] This implementation can not handle sparse chunks, try manually densifying them.
    .. [#densifies] This implementation densifies sparse chunks and therefore has increased memory usage.

    .. array-support:: pp.pca

    Parameters
    ----------
    data
        The (annotated) data matrix of shape `n_obs` × `n_vars`.
        Rows correspond to cells and columns to genes.
    n_comps
        Number of principal components to compute. Defaults to 50,
        or 1 - minimum dimension size of selected representation.
    layer
        If provided, which element of :attr:`~anndata.AnnData.layers` to use for PCA instead of `X`.
    obsm
        If provided, which element of :attr:`~anndata.AnnData.obsm` to use for PCA instead of `X`.
    zero_center
        If `True`, compute (or approximate) PCA from covariance matrix.
        If `False`, performa a truncated SVD instead of PCA.

        Our default PCA algorithms (see `svd_solver`) support implicit zero-centering,
        and therefore efficiently operating on sparse data.
    svd_solver
        SVD solver to use.
        See table above to see which solver class is used based on `chunked` and `zero_center`,
        as well as the default solver for each class when `svd_solver=None`.

        Efficient computation of the principal components of a sparse matrix
        currently only works with the `'arpack`' or `'covariance_eigh`' solver.

        `None`
            Choose automatically based on solver class (see table above).
        `'arpack'`
            ARPACK wrapper in SciPy (:func:`~scipy.sparse.linalg.svds`).
            Not available for *dask* arrays.
        `'covariance_eigh'`
            Classic eigendecomposition of the covariance matrix, suited for tall-and-skinny matrices.
            With dask, array must be CSR or dense and chunked as `(N, adata.shape[1])`.
        `'randomized'`
            Randomized algorithm from :cite:t:`Halko2009`.
            For *dask* arrays, this will use :func:`~dask.array.linalg.svd_compressed`.
        `'auto'`
            Choose automatically depending on the size of the problem:
            Will use `'full'` for small shapes and `'randomized'` for large shapes.
        `'tsqr'`
            “tall-and-skinny QR” algorithm from :cite:t:`Benson2013`.
            Only available for dense *dask* arrays.

        .. versionchanged:: 1.9.3
           Default value changed from `'arpack'` to None.
        .. versionchanged:: 1.4.5
           Default value changed from `'auto'` to `'arpack'`.
    chunked
        If `True`, perform an incremental PCA on segments of `chunk_size`.
        Automatically zero centers and ignores settings of `zero_center`, `random_seed` and `svd_solver`.
        If `False`, perform a full PCA/truncated SVD (see `svd_solver` and `zero_center`).
        See table above for which solver class is used.
    chunk_size
        Number of observations to include in each chunk.
        Required if `chunked=True` was passed.
    rng
        Change to use different initial states for the optimization.
    return_info
        Only relevant when not passing an :class:`~anndata.AnnData`:
        see “Returns”.
    {mask_var_hvg}
    layer
        Layer of `adata` to use as expression values.
    dtype
        Numpy data type string to which to convert the result.
    key_added
        If not specified, the embedding is stored as
        :attr:`~anndata.AnnData.obsm`\ `['X_pca']`, the loadings as
        :attr:`~anndata.AnnData.varm`\ `['PCs']`, and the the parameters in
        :attr:`~anndata.AnnData.uns`\ `['pca']`.
        If specified, the embedding is stored as
        :attr:`~anndata.AnnData.obsm`\ ``[key_added]``, the loadings as
        :attr:`~anndata.AnnData.varm`\ ``[key_added]``, and the the parameters in
        :attr:`~anndata.AnnData.uns`\ ``[key_added]``.
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

    `.obsm['X_pca' | key_added]` : :class:`~scipy.sparse.csr_matrix` | :class:`~scipy.sparse.csc_matrix` | :class:`~numpy.ndarray` (shape `(adata.n_obs, n_comps)`)
        PCA representation of data.
    `.varm['PCs' | key_added]` : :class:`~numpy.ndarray` (shape `(adata.n_vars, n_comps)`)
        The principal components containing the loadings *when `obsm=None`*.
    `.uns['pca' | key_added]['components']` : :class:`~numpy.ndarray` (shape `(adata.obsm[obsm].shape[1], n_comps)`)
        The principal components containing the loadings *when `obsm="..."`*.
    `.uns['pca' | key_added]['variance_ratio']` : :class:`~numpy.ndarray` (shape `(n_comps,)`)
        Ratio of explained variance.
    `.uns['pca' | key_added]['variance']` : :class:`~numpy.ndarray` (shape `(n_comps,)`)
        Explained variance, equivalent to the eigenvalues of the
        covariance matrix.

    """
    logg_start = logg.info("computing PCA")
    if (layer is not None or obsm is not None) and chunked:
        # Current chunking implementation relies on pca being called on X
        msg = "Cannot use `layer`/`obsm` and `chunked` at the same time."
        raise NotImplementedError(msg)

    # chunked calculation is not randomized, anyways
    if svd_solver in {"auto", "randomized"} and not chunked:
        logg.info(
            "Note that scikit-learn's randomized PCA might not be exactly "
            "reproducible across different computational platforms. For exact "
            "reproducibility, choose `svd_solver='arpack'`."
        )
    if return_anndata := isinstance(data, AnnData):
        if (layer is None and obsm is None) and not chunked and is_backed_type(data.X):
            msg = f"PCA is not implemented for matrices of type {type(data.X)} with chunked as False"
            raise NotImplementedError(msg)
        adata = data.copy() if copy else data
    else:
        adata = AnnData(data)

    # Unify new mask argument and deprecated use_highly_varible argument
    mask_var_param, mask_var = _handle_mask_var(
        adata, mask_var, obsm=obsm, use_highly_variable=use_highly_variable
    )
    del use_highly_variable
    adata_comp = adata[:, mask_var] if mask_var is not None else adata

    if n_comps is None:
        min_dim = min(adata_comp.n_vars, adata_comp.n_obs)
        n_comps = min_dim - 1 if min_dim <= settings.N_PCS else settings.N_PCS

    logg.info(f"    with {n_comps=}")

    x = _get_obs_rep(adata_comp, layer=layer, obsm=obsm)
    if is_backed_type(x) and (layer is not None or obsm is not None):
        msg = f"PCA is not implemented for matrices of type {type(x)} from layers/obsm"
        raise NotImplementedError(msg)

    # dask needs an int for random state
    if not isinstance(x, DaskArray):
        rng = np.random.default_rng(rng)
    elif not isinstance(rng, int):
        msg = f"rng needs to be an int, not a {type(rng).__name__} when passing a dask array"
        raise TypeError(msg)

    if chunked:
        if (
            not zero_center
            or rng is not None
            or (svd_solver is not None and svd_solver != "arpack")
        ):
            logg.debug("Ignoring zero_center, random_state, svd_solver")

        incremental_pca_kwargs = dict()
        if isinstance(x, DaskArray):
            from dask.array import zeros
            from dask_ml.decomposition import IncrementalPCA

            incremental_pca_kwargs["svd_solver"] = _handle_dask_ml_args(
                svd_solver, IncrementalPCA
            )
        else:
            from numpy import zeros
            from sklearn.decomposition import IncrementalPCA

        x_pca = zeros((x.shape[0], n_comps), x.dtype)

        pca_ = IncrementalPCA(n_components=n_comps, **incremental_pca_kwargs)

        for chunk, _, _ in adata_comp.chunked_X(chunk_size):
            chunk_dense = chunk.toarray() if isinstance(chunk, CSBase) else chunk
            pca_.partial_fit(chunk_dense)

        for chunk, start, end in adata_comp.chunked_X(chunk_size):
            chunk_dense = chunk.toarray() if isinstance(chunk, CSBase) else chunk
            x_pca[start:end] = pca_.transform(chunk_dense)
    elif zero_center:
        if isinstance(x, CSBase) and svd_solver == "lobpcg":
            msg = (
                f"{svd_solver=} for sparse relies on legacy code and will not be supported in the future. "
                "Also the lobpcg solver has been observed to be inaccurate. Please use 'arpack' instead."
            )
            warn(msg, FutureWarning)
            x_pca, pca_ = _pca_compat_sparse(x, n_comps, solver=svd_solver, rng=rng)
        else:
            if not isinstance(x, DaskArray):
                from sklearn.decomposition import PCA

                svd_solver = _handle_sklearn_args(
                    svd_solver, PCA, sparse=isinstance(x, CSBase)
                )
                pca_ = PCA(
                    n_components=n_comps,
                    svd_solver=svd_solver,
                    random_state=legacy_random_state(rng),
                )
            elif isinstance(x._meta, CSBase) or svd_solver == "covariance_eigh":
                from ._dask import PCAEighDask

                if rng is not None:
                    msg = f"Ignoring {rng=} when using a sparse dask array"
                    warn(msg, UserWarning)
                if svd_solver not in {None, "covariance_eigh"}:
                    msg = f"Ignoring {svd_solver=} when using a sparse dask array"
                    warn(msg, UserWarning)
                pca_ = PCAEighDask(n_components=n_comps)
            else:
                from dask_ml.decomposition import PCA

                svd_solver = _handle_dask_ml_args(svd_solver, PCA)
                pca_ = PCA(
                    n_components=n_comps,
                    svd_solver=svd_solver,
                    random_state=legacy_random_state(rng),
                )
            x_pca = pca_.fit_transform(x)
    else:
        if isinstance(x, DaskArray):
            if isinstance(x._meta, CSBase):
                msg = (
                    "`zero_center=False` is not supported for sparse Dask arrays (yet). "
                    "See <https://github.com/dask/dask-ml/issues/123>."
                )
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
            n_components=n_comps,
            random_state=legacy_random_state(rng),
            algorithm=svd_solver,
        )
        x_pca = pca_.fit_transform(x)

    if x_pca.dtype.descr != np.dtype(dtype).descr:
        x_pca = x_pca.astype(dtype)

    if return_anndata:
        key_obsm, key_varm, key_uns = (
            ("X_pca", "PCs", "pca") if key_added is None else [key_added] * 3
        )
        adata.obsm[key_obsm] = x_pca

        if obsm:
            pass  # see below, components are stored in `uns`.
        elif mask_var is not None:
            adata.varm[key_varm] = np.zeros(shape=(adata.n_vars, n_comps))
            adata.varm[key_varm][mask_var] = pca_.components_.T
        else:
            adata.varm[key_varm] = pca_.components_.T

        adata.uns[key_uns] = dict(
            params=dict(
                zero_center=zero_center,
                use_highly_variable=mask_var_param == "highly_variable",
                mask_var=mask_var_param,
                **(dict(layer=layer) if layer is not None else {}),
                **(dict(obsm=obsm) if obsm is not None else {}),
            ),
            variance=pca_.explained_variance_,
            variance_ratio=pca_.explained_variance_ratio_,
            **(dict(components=pca_.components_.T) if obsm is not None else {}),
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
                x_pca,
                pca_.components_,
                pca_.explained_variance_ratio_,
                pca_.explained_variance_,
            )
        else:
            return x_pca


def _handle_mask_var(
    adata: AnnData,
    mask_var: NDArray[np.bool_] | str | Empty | None,
    *,
    obsm: str | None = None,
    use_highly_variable: bool | None,
) -> tuple[np.ndarray | str | None, np.ndarray | None]:
    """Unify new mask argument and deprecated use_highly_varible argument.

    Returns both the normalized mask parameter and the validated mask array.
    """
    if obsm:
        if mask_var is not _empty and mask_var is not None:
            msg = "Argument `mask_var` is incompatible with `obsm`."
            raise ValueError(msg)
        return None, None

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


def _handle_x_args[T: LiteralString](
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
        warn(msg, UserWarning)
    return default
