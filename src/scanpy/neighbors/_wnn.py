"""Weighted nearest neighbors (WNN) across modalities."""

from __future__ import annotations

import itertools
import operator
from collections.abc import Mapping, Sequence
from typing import TYPE_CHECKING

import numba
import numpy as np
from fast_array_utils.numba import njit
from scipy import sparse

from .. import logging as logg
from .._compat import CSBase, warn
from .._docs import doc_rng
from .._utils import _doc_params, get_literal_vals
from .._utils._doctests import doctest_needs
from .._utils.random import _accepts_legacy_random_state, _LegacyRng
from ..get.get import _rep_to_json
from . import _connectivity
from ._common import (
    _get_indices_distances_from_rect_matrix,
    _get_sparse_matrix_from_indices_distances,
    _make_transformer,
)
from ._types import _KnownTransformer, _Method

if TYPE_CHECKING:
    from typing import Any

    from anndata import AnnData
    from numpy.typing import NDArray

    from .._compat import CSRBase
    from .._utils.random import RNGLike, SeedLike
    from ..get.get import RepAcc
    from ._types import KnnTransformerLike


@doctest_needs("igraph")
@_accepts_legacy_random_state(0)
@_doc_params(rng=doc_rng)
def wnn(  # noqa: PLR0912, PLR0913, PLR0915
    adata: AnnData,
    use_rep: Sequence[RepAcc | str] | Mapping[str, RepAcc | str],
    n_pcs: int | Sequence[int | None] | None = None,
    *,
    n_neighbors: int = 20,
    knn_range: int = 200,
    l2_norm: bool = True,
    sd_scale: float = 1.0,
    cross_constant: float | Sequence[float] = 1e-4,
    smooth: bool = False,
    prune_snn: float = 1 / 15,
    transformer: KnnTransformerLike | _KnownTransformer | None = None,
    method: _Method = "umap",
    rng: SeedLike | RNGLike | None = None,
    key_added: str | None = None,
    weight_key: str | None = "{}_weight",
    copy: bool = False,
) -> AnnData | None:
    r"""Compute a weighted nearest neighbor (WNN) graph across modalities :cite:p:`Hao2021`.

    Implementation of Seurat’s ``FindMultiModalNeighbors``.
    For every cell, each modality gets a weight from how well the cell’s neighbors
    in that modality predict its embedding, compared to the neighbors from the other modalities.
    The neighbors are then chosen by a weighted sum of per-modality affinity kernels,
    whose widths are set from the shared nearest neighbor structure, as in Seurat.

    Use this as an alternative to :func:`~scanpy.pp.neighbors`:
    it writes the same fields, so downstream steps
    (e.g. :func:`~scanpy.tl.umap` or :func:`~scanpy.tl.leiden`) work unchanged.
    This CPU implementation is based on the rapids-singlecell package.

    .. array-support:: pp.wnn

    Parameters
    ----------
    adata
        Annotated data matrix, or a :class:`~mudata.MuData` object
        whose modalities share the same cells.
    use_rep
        One dimensional reduction per modality, e.g. ``['X_pca', 'X_apca']``
        (keys of `.obsm` or accessors, see :func:`~scanpy.pp.neighbors`).
        For a :class:`~mudata.MuData` object, a mapping of modality to a representation
        of that modality, e.g. ``{{'rna': 'X_pca', 'prot': 'X_pca'}}``.
        A mapping also names the modalities for `weight_key` with :class:`~anndata.AnnData`.
    n_pcs
        Number of dimensions to use, per modality or for all of them
        (Seurat’s ``dims.list``). `None` uses all dimensions.
    n_neighbors
        Number of multimodal neighbors (Seurat’s ``k.nn``),
        also used to compute the modality weights.
        Unlike :func:`~scanpy.pp.neighbors`, this excludes the cell itself,
        so ``.uns[key_added]['params']['n_neighbors']`` is `n_neighbors` + 1.
    knn_range
        Number of neighbors searched per modality,
        whose union are the candidates for the multimodal neighbors (Seurat’s ``knn.range``).
    l2_norm
        L2-normalize the embeddings of every cell.
    sd_scale
        Scaling factor of the kernel widths.
    cross_constant
        Constant to avoid dividing by zero when comparing the within- and
        cross-modality predictions, per modality or for all of them
        (Seurat’s ``cross.contant.list``).
    smooth
        Average the modality scores over each cell’s neighbors in that modality.
        Seurat’s ``smooth=TRUE`` fails (as of Seurat 5.5.1);
        this follows its intent, averaging every pair of modalities’ scores.
    prune_snn
        Jaccard cutoff for edges of the shared nearest neighbor graph (Seurat’s ``prune.SNN``).
    transformer
        kNN search backend used within each modality, following the API of
        :class:`~sklearn.neighbors.KNeighborsTransformer`.
        It needs to use the euclidean metric
        and return at least `knn_range` neighbors (including the cell itself).
        See :doc:`/how-to/knn-transformers` for more details.
        Also accepts the following known options:

        `None` (the default)
            Behavior depends on data size.
            For small data, we will calculate exact kNN, otherwise we use
            :class:`~pynndescent.pynndescent_.PyNNDescentTransformer`
            (Seurat uses an approximate Annoy search).
        `'sklearn'`
            Exact :class:`~sklearn.neighbors.KNeighborsTransformer`.
        `'pynndescent'`
            :class:`~pynndescent.pynndescent_.PyNNDescentTransformer`
    method
        Kernel that derives `connectivities` from the multimodal neighbors,
        see :func:`~scanpy.pp.neighbors`.
    {rng}

        *ignored if* ``transformer`` *is an instance.*
    key_added
        If not specified, the neighbors data is stored in `.uns['neighbors']`,
        distances, connectivities and the SNN graph in
        `.obsp['distances']`, `.obsp['connectivities']` and `.obsp['snn']`.
        If specified, the neighbors data is added to `.uns[key_added]`,
        and the graphs to `.obsp[f'{{key_added}}_distances']`,
        `.obsp[f'{{key_added}}_connectivities']` and `.obsp[f'{{key_added}}_snn']`.
    weight_key
        Format string for the `.obs` columns holding the modality weights,
        filled with the modality name (the `use_rep` string or mapping key).
        ``'{{}}:mod_weight'`` gives muon’s column names, `None` skips the weights.
    copy
        Return a copy instead of writing to `adata`.

    Returns
    -------
    Returns `None` if `copy=False`, else returns an `AnnData` object. Sets the following fields:

    `adata.obsp['distances' | f'{{key_added}}_distances']` : :class:`scipy.sparse.csr_matrix` (dtype `float`)
        The `n_neighbors` multimodal neighbors of every cell, excluding the cell itself,
        with Seurat’s WNN distance :math:`\sqrt{{(1 - s)/2}}` of the weighted kernel similarity :math:`s`.
        Rows keep the neighbors in decreasing similarity, as Seurat does.
    `adata.obsp['connectivities' | f'{{key_added}}_connectivities']` : :class:`scipy.sparse.csr_matrix` (dtype `float`)
        Weighted adjacency matrix of the multimodal neighbors computed with `method`,
        counting the cell itself as its first neighbor (i.e. ``n_neighbors + 1``).
    `adata.obsp['snn' | f'{{key_added}}_snn']` : :class:`scipy.sparse.csr_matrix` (dtype `float32`)
        Seurat’s ``wsnn`` graph, the Jaccard index of the multimodal neighbor sets,
        with Seurat’s unit diagonal. Seurat clusters on this graph,
        e.g. ``sc.tl.leiden(adata, obsp='snn')``.
    `adata.obs[weight_key.format(name)]` : :class:`~pandas.Series` (dtype `float`)
        The weight of each modality per cell.
    `adata.uns['neighbors' | key_added]` : :class:`dict`
        neighbors parameters.
        The representation and number of dimensions of each modality are stored in
        ``params['modalities']``.
        For :class:`~anndata.AnnData`, ``params['use_rep']`` and ``params['n_pcs']``
        refer to the first modality (used by :func:`~scanpy.tl.umap` for initialization).

    Examples
    --------
    >>> import numpy as np
    >>> import scanpy as sc
    >>> adata = sc.datasets.pbmc68k_reduced()
    >>> rng = np.random.default_rng(0)
    >>> adata.obsm["X_apca"] = adata.obsm["X_pca"][:, :10] + rng.normal(
    ...     size=(adata.n_obs, 10)
    ... )
    >>> sc.pp.wnn(adata, {{"RNA": "X_pca", "ADT": "X_apca"}}, knn_range=100)
    >>> sc.tl.leiden(adata, obsp="snn", flavor="igraph")
    >>> sc.tl.umap(adata)

    See Also
    --------
    :func:`~scanpy.pp.neighbors`

    """
    start = logg.info("computing weighted nearest neighbors")

    # validate before any work, so that errors leave `adata` unchanged
    if method not in (methods := get_literal_vals(_Method)):
        msg = f"`method` needs to be one of {methods}."
        raise ValueError(msg)
    _check_transformer(transformer)
    is_mudata = hasattr(adata, "mod")
    names, use_reps, n_pcs = _parse_use_rep(adata, use_rep, n_pcs, is_mudata=is_mudata)
    weight_cols = [] if weight_key is None else [weight_key.format(n) for n in names]
    if len(set(weight_cols)) < len(weight_cols) or "" in weight_cols:
        msg = (
            f"`weight_key` {weight_key!r} must give a distinct column per modality, "
            "e.g. '{}_weight'."
        )
        raise ValueError(msg)
    n_mod = len(names)
    n_obs = adata.n_obs
    n_neighbors, knn_range = operator.index(n_neighbors), operator.index(knn_range)
    if n_neighbors < 2:
        msg = "`n_neighbors` must be at least 2."
        raise ValueError(msg)
    if knn_range <= n_neighbors:
        msg = "`knn_range` must be larger than `n_neighbors`."
        raise ValueError(msg)
    if knn_range > n_obs:
        msg = f"`knn_range` must be at most the number of cells ({n_obs})."
        raise ValueError(msg)
    if not (np.isfinite(sd_scale) and sd_scale > 0):
        msg = "`sd_scale` must be positive."
        raise ValueError(msg)
    if not 0 <= prune_snn <= 1:
        msg = "`prune_snn` must be between 0 and 1."
        raise ValueError(msg)
    cross_constant = np.asarray(cross_constant, dtype=np.float64)
    if cross_constant.ndim == 0:
        cross_constant = np.full(n_mod, cross_constant)
    elif cross_constant.shape != (n_mod,):
        msg = "`cross_constant` needs one entry per modality."
        raise ValueError(msg)
    if not (np.isfinite(cross_constant).all() and (cross_constant > 0).all()):
        msg = "`cross_constant` must be positive."
        raise ValueError(msg)
    meta_random_state = (
        dict(random_state=rng.arg) if isinstance(rng, _LegacyRng) else {}
    )
    rng = np.random.default_rng(rng)
    # the default representation (`None`) may run PCA, which belongs on the copy
    copied = copy and any(rep is None for rep in use_reps)
    if copied:
        adata = adata.copy()
    modalities, reps = _get_reps(adata, names, use_reps, n_pcs, is_mudata=is_mudata)

    if l2_norm:
        reps = [_l2_normalize(x) for x in reps]

    # kNN within each modality: the first `n_neighbors` give the weights,
    # all `knn_range` the candidates for the multimodal neighbors.
    knns = np.empty((n_mod, n_obs, knn_range), dtype=np.int32)
    for m, x in enumerate(reps):
        knns[m] = _knn(
            x, knn_range=knn_range, transformer=transformer, rng=rng, name=names[m]
        )
    logg.debug("    computed neighbors within each modality")

    # dist[m, r]: distance of modality m to its prediction from modality r’s neighbors
    # the same distance kernel as `_multimodal_knn`, so that a candidate that is
    # every modality’s nearest neighbor gets a WNN distance of exactly 0
    emb = np.concatenate(reps, axis=1)
    dim_off = np.cumsum([0] + [x.shape[1] for x in reps])
    nearest = _nearest_dist(emb, dim_off, knns)
    sigma = np.empty((n_mod, n_obs))
    dist = np.empty((n_mod, n_mod, n_obs))
    for m, x in enumerate(reps):
        for r in range(n_mod):
            dist[m, r] = _impute_dist(x, knns[r], n_neighbors, nearest[m])
        post_off, postings = _postings(knns[m], n_neighbors)
        sigma[m] = _snn_bandwidth(
            knns[m], n_neighbors, post_off, postings, x, nearest[m]
        )
    if (n_zero := (sigma == 0).sum(axis=1)).any():
        counts = ", ".join(
            f"{c} in {n!r}" for n, c in zip(names, n_zero, strict=True) if c
        )
        msg = (
            f"Cells ({counts}) have the same embedding as all of their nearest "
            "neighbors in that modality, so its kernel width is zero "
            "and their modality weights are unreliable (Seurat gives NaN). "
            "Remove duplicated or empty cells before running WNN."
        )
        warn(msg, UserWarning)
    sigma *= sd_scale
    # keep the weights finite for zero widths
    sigma = np.maximum(sigma, np.finfo(np.float64).tiny)

    kernel = np.exp(-dist / sigma[:, None, :])
    exp_scores = np.zeros((n_mod, n_obs))
    for m, r in itertools.permutations(range(n_mod), 2):
        score = np.clip(kernel[m, m] / (kernel[m, r] + cross_constant[m]), 0, 200)
        if smooth:
            score = score[knns[m, :, 1:n_neighbors]].mean(axis=1)
        exp_scores[m] += np.exp(score)
    weights = exp_scores / exp_scores.sum(axis=0)
    logg.debug("    computed modality weights")

    knn_indices, knn_dissim = _multimodal_knn(
        emb, dim_off, knns, n_neighbors, weights, sigma, nearest
    )
    # Seurat’s sqrt((1 - s) / 2) of the weighted kernel similarity s
    knn_dist = np.sqrt(np.clip(knn_dissim / 2, 0, 1))
    del emb, knns

    # the connectivity methods expect the cell itself in column 0
    knn_indices = np.column_stack((np.arange(n_obs), knn_indices))
    knn_dist = np.column_stack((np.zeros(n_obs), knn_dist))
    distances = _get_sparse_matrix_from_indices_distances(
        knn_indices, knn_dist, keep_self=False
    )
    post_off, postings = _postings(knn_indices[:, 1:], n_neighbors)
    snn = _snn_graph(knn_indices[:, 1:], post_off, postings, prune=prune_snn)
    connectivities = _connectivities(
        method, knn_indices, knn_dist, distances, n_neighbors=n_neighbors + 1
    )

    if copy and not copied:
        adata = adata.copy()
    if adata.is_view:
        adata._init_as_actual(adata.copy())
    # `tl.umap` uses one representation to initialize disconnected components,
    # MuData’s are per modality, where it falls back to the default.
    layout_rep = {} if is_mudata else modalities[names[0]]
    key_added = "neighbors" if key_added is None else key_added
    prefix = "" if key_added == "neighbors" else f"{key_added}_"
    adata.uns[key_added] = dict(
        connectivities_key=f"{prefix}connectivities",
        distances_key=f"{prefix}distances",
        snn_key=f"{prefix}snn",
        params=dict(
            # scanpy’s convention counts the cell itself
            n_neighbors=n_neighbors + 1,
            method=method,
            **meta_random_state,
            metric="euclidean",
            **layout_rep,
            modalities=modalities,
            knn_range=knn_range,
            l2_norm=l2_norm,
            sd_scale=sd_scale,
            cross_constant=[float(c) for c in cross_constant],
            smooth=smooth,
            prune_snn=prune_snn,
        ),
    )
    adata.obsp[f"{prefix}distances"] = distances
    adata.obsp[f"{prefix}connectivities"] = connectivities
    adata.obsp[f"{prefix}snn"] = snn
    if weight_key is not None:
        for col, w in zip(weight_cols, weights, strict=True):
            adata.obs[col] = w

    logg.info(
        "    finished",
        time=start,
        deep=(
            f"added to `.uns[{key_added!r}]`\n"
            f"    `.obsp[{prefix + 'distances'!r}]`, distances for each pair of neighbors\n"
            f"    `.obsp[{prefix + 'connectivities'!r}]`, weighted adjacency matrix\n"
            f"    `.obsp[{prefix + 'snn'!r}]`, shared nearest neighbor graph"
            + "".join(f"\n    `.obs[{col!r}]`, modality weights" for col in weight_cols)
        ),
    )
    return adata if copy else None


def _parse_use_rep(
    adata: AnnData,
    use_rep: Sequence[RepAcc | str] | Mapping[str, RepAcc | str],
    n_pcs: int | Sequence[int | None] | None,
    *,
    is_mudata: bool,
) -> tuple[list[str], list[RepAcc | str | None], list[int | None]]:
    """Get the modality names, representations, and numbers of dimensions."""
    if isinstance(use_rep, Mapping):
        names, reps = list(use_rep), list(use_rep.values())
    elif is_mudata:
        msg = (
            "For MuData, pass `use_rep` as a mapping of modality to representation, "
            "e.g. `{'rna': 'X_pca', 'prot': 'X_pca'}`."
        )
        raise TypeError(msg)
    # accessors are iterable without being sequences
    elif isinstance(use_rep, Sequence | np.ndarray) and not isinstance(use_rep, str):
        names = reps = list(use_rep)
    else:
        msg = "`use_rep` needs one representation per modality, e.g. `['X_pca', 'X_apca']`."
        raise TypeError(msg)
    if len(reps) < 2:
        msg = "WNN needs at least two modalities."
        raise ValueError(msg)
    if len(set(names)) != len(names) or not all(isinstance(n, str) for n in names):
        msg = (
            f"Modality names {names} must be unique strings; pass `use_rep` as a "
            "mapping, e.g. `{'pca10': 'X_pca', 'pca30': 'X_pca'}`."
        )
        raise ValueError(msg)
    if n_pcs is None or isinstance(n_pcs, int | np.integer):
        n_pcs = [n_pcs] * len(reps)
    elif len(n_pcs) != len(reps):
        msg = "`n_pcs` needs one entry per modality."
        raise ValueError(msg)
    if is_mudata:
        for name in names:
            if not adata.mod[name].obs_names.equals(adata.obs_names):
                msg = (
                    f"Modality {name!r} does not hold the cells of the MuData object "
                    "in the same order. Run `muon.pp.intersect_obs` and reorder it, "
                    f"e.g. `mdata.mod[{name!r}] = mdata[{name!r}][mdata.obs_names].copy()`."
                )
                raise ValueError(msg)
    return names, reps, list(n_pcs)


def _get_reps(
    adata: AnnData,
    names: Sequence[str],
    use_reps: Sequence[RepAcc | str | None],
    n_pcs: Sequence[int | None],
    *,
    is_mudata: bool,
) -> tuple[dict[str, dict[str, Any]], list[NDArray[np.float64]]]:
    """Get the per-modality metadata and representations."""
    from ..tools._utils import _choose_representation_compat

    modalities, xs = {}, []
    for name, rep, npc in zip(names, use_reps, n_pcs, strict=True):
        mod = adata.mod[name] if is_mudata else adata
        x = _choose_representation_compat(mod, use_rep=rep, n_pcs=npc)
        if isinstance(x, CSBase):
            msg = f"Representation {rep!r} must be dense."
            raise TypeError(msg)
        xs.append(np.ascontiguousarray(x, dtype=np.float64))
        modalities[name] = dict(use_rep=_rep_to_json(rep), n_pcs=x.shape[1])
    return modalities, xs


def _check_transformer(
    transformer: KnnTransformerLike | _KnownTransformer | None,
) -> None:
    if transformer is None:
        return
    if isinstance(transformer, str):
        if transformer not in (known := get_literal_vals(_KnownTransformer)):
            msg = f"`transformer` needs to be one of {known} or an instance."
            raise ValueError(msg)
        return
    # WNN’s kernels are euclidean, so the neighbors need to be, too
    metric = getattr(transformer, "metric", "euclidean")
    if metric == "minkowski" and getattr(transformer, "p", 2) == 2:
        return
    if metric not in {"euclidean", "l2"}:
        msg = f"`transformer` needs to use the euclidean metric, not {metric!r}."
        raise ValueError(msg)


def _l2_normalize(x: NDArray[np.float64]) -> NDArray[np.float64]:
    norm = np.linalg.norm(x, axis=1, keepdims=True)
    # as Seurat’s `L2Norm`, all-zero rows stay zero
    return x / np.where(norm > 0, norm, 1)


def _knn(
    x: NDArray[np.float64],
    *,
    knn_range: int,
    transformer: KnnTransformerLike | _KnownTransformer | None,
    rng: np.random.Generator,
    name: str,
) -> NDArray[np.int32]:
    """Find the `knn_range` nearest neighbors of every cell, the cell itself first."""
    from sklearn.base import clone

    n_obs = x.shape[0]
    shortcut = transformer == "sklearn" or (transformer is None and n_obs < 8192)
    proto = _make_transformer(
        transformer,
        shortcut=shortcut,
        n_index=n_obs,
        n_neighbors=knn_range,
        metric="euclidean",
        metric_params={},
        rng=rng,
    )
    knn = clone(proto)
    if shortcut:
        # ask for exactly `knn_range`: `transform` would add one
        d = knn.fit(x).kneighbors_graph(x, n_neighbors=knn_range, mode="distance")
    else:
        d = knn.fit_transform(x)
    try:
        indices, distances = _get_indices_distances_from_rect_matrix(d, knn_range)
    except ValueError as e:
        msg = (
            f"The kNN search for modality {name!r} returned too few neighbors: {e} "
            "Use a `transformer` that returns at least `knn_range` neighbors."
        )
        raise ValueError(msg) from None
    # the numba kernels index with these unchecked
    if indices.min() < 0 or indices.max() >= n_obs or not np.isfinite(distances).all():
        msg = (
            f"The kNN search for modality {name!r} returned invalid neighbors "
            "(indices out of range or non-finite distances)."
        )
        raise ValueError(msg)
    return _ensure_self_first(indices.astype(np.int32))


def _ensure_self_first(knn: NDArray[np.int32]) -> NDArray[np.int32]:
    """Put every cell first in its own neighbor list.

    Seurat assumes column 0 is the cell itself.
    Exact duplicates can push it further back (or out of the list),
    so move it to the front, keeping the order of the remaining neighbors.
    """
    n, k = knn.shape
    rows = np.arange(n, dtype=knn.dtype)
    bad = np.flatnonzero(knn[:, 0] != rows)
    if bad.size == 0:
        return knn
    sub = knn[bad]
    is_self = sub == rows[bad, None]
    has_self = is_self.any(axis=1)
    order = np.argsort(np.where(is_self, -1, np.arange(k)), axis=1, kind="stable")
    moved = np.take_along_axis(sub, order, axis=1)
    shifted = np.concatenate((rows[bad, None], sub[:, :-1]), axis=1)
    knn[bad] = np.where(has_self[:, None], moved, shifted)
    return knn


@njit
def _nearest_dist(
    emb: NDArray[np.float64], dim_off: NDArray[np.int64], knns: NDArray[np.int32]
) -> NDArray[np.float64]:
    """Distance of every cell to its nearest neighbor, per modality."""
    n_mod, n, _ = knns.shape
    out = np.empty((n_mod, n))
    for i in numba.prange(n):
        for m in range(n_mod):
            xi = emb[i, dim_off[m] : dim_off[m + 1]]
            xc = emb[knns[m, i, 1], dim_off[m] : dim_off[m + 1]]
            out[m, i] = np.sqrt(np.sum((xi - xc) ** 2))
    return out


@njit
def _impute_dist(
    emb: NDArray[np.float64],
    knn: NDArray[np.int32],
    k: int,
    nearest: NDArray[np.float64],
) -> NDArray[np.float64]:
    """Seurat’s ``PredictAssay`` + ``impute_dist``.

    Distance between each cell’s embedding and the mean embedding of its neighbors
    ``1..k-1`` (dropping the cell itself), minus the distance to the nearest neighbor,
    clipped at zero.
    """
    n, d = emb.shape
    out = np.empty(n)
    for i in numba.prange(n):
        pred = np.zeros(d)
        for c in range(1, k):
            pred += emb[knn[i, c]]
        pred /= k - 1
        out[i] = max(np.sqrt(np.sum((emb[i] - pred) ** 2)) - nearest[i], 0.0)
    return out


def _postings(
    knn: NDArray[np.int32 | np.int64], s: int
) -> tuple[NDArray[np.int64], NDArray[np.int32]]:
    """Reverse neighbor lists of the first `s` neighbors.

    Cell `j`’s list ``postings[post_off[j]:post_off[j + 1]]``
    contains the cells that have `j` among their first `s` neighbors.
    """
    flat = knn[:, :s].ravel()
    post_off = np.zeros(knn.shape[0] + 1, dtype=np.int64)
    np.cumsum(np.bincount(flat, minlength=knn.shape[0]), out=post_off[1:])
    postings = (np.argsort(flat, kind="stable") // s).astype(np.int32)
    return post_off, postings


@numba.njit(cache=True)  # noqa: TID251
def _snn_partners(
    knn: NDArray[np.int32 | np.int64],
    s: int,
    i: int,
    post_off: NDArray[np.int64],
    postings: NDArray[np.int32],
) -> tuple[NDArray[np.int32], NDArray[np.int64]]:
    """Shared nearest neighbors of cell `i`, in ascending order, and their overlap.

    The partners of `i` are the reverse neighbor lists of its first `s` neighbors:
    a cell sharing `c` neighbors with `i` appears `c` times,
    the overlap Seurat’s ``ComputeSNN`` counts.
    """
    n_keys = 0
    for t in range(s):
        n_keys += post_off[knn[i, t] + 1] - post_off[knn[i, t]]
    keys = np.empty(n_keys, dtype=np.int32)
    at = 0
    for t in range(s):
        beg, end = post_off[knn[i, t]], post_off[knn[i, t] + 1]
        keys[at : at + end - beg] = postings[beg:end]
        at += end - beg
    keys.sort()
    partners = np.empty(n_keys, dtype=np.int32)
    counts = np.zeros(n_keys, dtype=np.int64)
    n_partners = 0
    for p in range(n_keys):
        if p == 0 or keys[p] != keys[p - 1]:
            partners[n_partners] = keys[p]
            n_partners += 1
        counts[n_partners - 1] += 1
    return partners[:n_partners], counts[:n_partners]


@njit
def _snn_bandwidth(  # noqa: PLR0917
    knn: NDArray[np.int32],
    s: int,
    post_off: NDArray[np.int64],
    postings: NDArray[np.int32],
    emb: NDArray[np.float64],
    nearest: NDArray[np.float64],
) -> NDArray[np.float64]:
    """Seurat’s ``ComputeSNNwidth`` (``SNN_SmallestNonzero_Dist``).

    The mean distance (minus the nearest neighbor distance, clipped at zero)
    to the `s` farthest of the cell’s least shared neighbors,
    i.e. those whose Jaccard index is at most the `s`-th smallest one.
    """
    n = knn.shape[0]
    sigma = np.empty(n)
    for i in numba.prange(n):
        partners, counts = _snn_partners(knn, s, i, post_off, postings)
        jaccard = counts / (2 * s - counts)
        n_i = min(s, partners.size)
        thresh = np.sort(jaccard)[n_i - 1]
        dist = np.empty(partners.size)
        n_sel = 0
        for p in range(partners.size):
            if jaccard[p] <= thresh:
                d = np.sqrt(np.sum((emb[i] - emb[partners[p]]) ** 2))
                dist[n_sel] = max(d - nearest[i], 0.0)
                n_sel += 1
        sigma[i] = np.sort(dist[:n_sel])[n_sel - n_i :].mean()
    return sigma


@njit
def _snn_counts(
    knn: NDArray[np.int64],
    post_off: NDArray[np.int64],
    postings: NDArray[np.int32],
    prune: float,
) -> NDArray[np.int64]:
    n, k = knn.shape
    out = np.empty(n, dtype=np.int64)
    for i in numba.prange(n):
        _, counts = _snn_partners(knn, k, i, post_off, postings)
        out[i] = np.sum(counts / (2 * k - counts) >= prune)
    return out


@njit
def _snn_fill(  # noqa: PLR0917
    knn: NDArray[np.int64],
    post_off: NDArray[np.int64],
    postings: NDArray[np.int32],
    prune: float,
    indptr: NDArray[np.int64],
    indices: NDArray[np.int32],
    data: NDArray[np.float32],
) -> None:
    n, k = knn.shape
    for i in numba.prange(n):
        partners, counts = _snn_partners(knn, k, i, post_off, postings)
        jaccard = counts / (2 * k - counts)
        at = indptr[i]
        for p in range(partners.size):
            if jaccard[p] >= prune:
                indices[at] = partners[p]
                data[at] = jaccard[p]
                at += 1


def _snn_graph(
    knn: NDArray[np.int64],
    post_off: NDArray[np.int64],
    postings: NDArray[np.int32],
    *,
    prune: float,
) -> CSRBase:
    """Seurat’s ``ComputeSNN``: Jaccard index of the neighbor sets, pruned."""
    n = knn.shape[0]
    indptr = np.zeros(n + 1, dtype=np.int64)
    np.cumsum(_snn_counts(knn, post_off, postings, prune), out=indptr[1:])
    indices = np.empty(indptr[-1], dtype=np.int32)
    data = np.empty(indptr[-1], dtype=np.float32)
    _snn_fill(knn, post_off, postings, prune, indptr, indices, data)
    if indptr[-1] <= np.iinfo(np.int32).max:
        indptr = indptr.astype(np.int32)
    snn = sparse.csr_matrix((data, indices, indptr), shape=(n, n))  # noqa: TID251
    snn.has_canonical_format = True
    return snn


@njit
def _multimodal_knn(  # noqa: PLR0917
    emb: NDArray[np.float64],
    dim_off: NDArray[np.int64],
    knns: NDArray[np.int32],
    k: int,
    weights: NDArray[np.float64],
    sigma: NDArray[np.float64],
    nearest: NDArray[np.float64],
) -> tuple[NDArray[np.int64], NDArray[np.float64]]:
    """Seurat’s ``MultiModalNN``.

    The candidates are the union of each modality’s neighbors (excluding the cell itself),
    in order of first occurrence.
    They are ranked by the weighted sum of the per-modality kernels,
    keeping the candidate order for ties, as R’s ``order`` does.
    Returns one minus that similarity, computed so that it is exactly 0
    for a candidate that is the nearest neighbor in every modality.
    """
    n_mod, n, _ = knns.shape
    out_idx = np.empty((n, k), dtype=np.int64)
    out_dissim = np.empty((n, k))
    for i in numba.prange(n):
        cand = knns[:, i, 1:].flatten()
        # keep the first occurrence of every candidate
        order = np.argsort(cand, kind="mergesort")
        keep = np.ones(cand.size, dtype=np.bool)
        for p in range(1, cand.size):
            if cand[order[p]] == cand[order[p - 1]]:
                keep[order[p]] = False
        cand = cand[keep]
        dissim = np.zeros(cand.size)
        for m in range(n_mod):
            xi = emb[i, dim_off[m] : dim_off[m + 1]]
            for c in range(cand.size):
                xc = emb[cand[c], dim_off[m] : dim_off[m + 1]]
                d = max(np.sqrt(np.sum((xi - xc) ** 2)) - nearest[m, i], 0.0)
                # the weights sum to 1, so this is 1 - sum(w * exp(-d / sigma))
                dissim[c] -= weights[m, i] * np.expm1(-d / sigma[m, i])
        best = np.argsort(dissim, kind="mergesort")[:k]
        out_idx[i] = cand[best]
        out_dissim[i] = dissim[best]
    return out_idx, out_dissim


def _connectivities(
    method: _Method,
    knn_indices: NDArray[np.int64],
    knn_dist: NDArray[np.float64],
    distances: CSRBase,
    *,
    n_neighbors: int,
) -> CSRBase:
    n_obs = knn_indices.shape[0]
    if method == "umap":
        return _connectivity.umap(
            knn_indices, knn_dist, n_obs=n_obs, n_neighbors=n_neighbors
        )
    if method == "gauss":
        # WNN distances depend on the direction (on the weights and kernel widths
        # of the cell whose neighbors they are), so average both directions
        # `gauss` sorts the indices of its input in place
        connectivities = _connectivity.gauss(distances.copy(), n_neighbors, knn=True)
        return ((connectivities + connectivities.T) / 2).tocsr()
    return _connectivity.jaccard(knn_indices, n_obs=n_obs, n_neighbors=n_neighbors)
