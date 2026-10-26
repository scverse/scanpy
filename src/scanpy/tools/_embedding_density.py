"""Calculate density of cells in embeddings."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
from scverse_misc import Deprecation, deprecated_arg

from .. import logging as logg
from .._docs import DEPR_KEY_ADDED, doc_out_key
from .._settings import Default
from .._utils import _doc_params, _get_basis_key, sanitize_anndata
from ..get._out import _out_key
from ..get._types import MultiAcc
from ..get._use import _obs_col, _resolve_ref

if TYPE_CHECKING:
    from collections.abc import Sequence

    from anndata import AnnData
    from anndata.acc import AdRef


def _calc_density(x: np.ndarray, y: np.ndarray):
    """Calculate the density of points in 2 dimensions."""
    from scipy.stats import gaussian_kde

    # Calculate the point density
    xy = np.vstack([x, y])
    z = gaussian_kde(xy)(xy)

    min_z = np.min(z)
    max_z = np.max(z)

    # Scale between 0 and 1
    scaled_z = (z - min_z) / (max_z - min_z)

    return scaled_z


@_doc_params(
    out=doc_out_key(
        "the densities",
        example="A.obs['density']",
        key_added="`key_added='k'` is `out=A.obs['k']`.",
        default="`.obs['[basis]_density_[groupby]']` and `.uns['[basis]_density_[groupby]_params']`",
    ).replace("are stored under the same key", "are named after its key")
)
@deprecated_arg("key_added", DEPR_KEY_ADDED)
@deprecated_arg("basis", Deprecation("1.13.0", "Use `use` instead."))
def embedding_density(  # noqa: PLR0912, PLR0915
    adata: AnnData,
    use: MultiAcc | str | Default = Default("'umap'"),
    *,
    groupby: AdRef | str | None = None,
    components: str | Sequence[str] | None = None,
    out: AdRef | str | bool = False,
    # deprecated
    basis: MultiAcc | str | None = None,
    key_added: str | None = None,
) -> np.ndarray | None:
    r"""Calculate the density of cells in an embedding (per condition).

    Gaussian kernel density estimation is used to calculate the density of
    cells in an embedded space. This can be performed per category over a
    categorical cell annotation. The cell density can be plotted using the
    `pl.embedding_density` function.

    Note that density values are scaled to be between 0 and 1. Thus, the
    density value at each cell is only comparable to densities in
    the same category.

    Beware that the KDE estimate used (`scipy.stats.gaussian_kde`) becomes
    unreliable if you don't have enough cells in a category.

    This function was written by Sophie Tritschler and implemented into
    Scanpy by Malte Luecken.

    .. array-support:: tl.embedding_density

    Parameters
    ----------
    adata
        The annotated data matrix.
    use
        The embedding over which the density will be calculated,
        e.g. `A.obsm['umap']`.
        If :attr:`scanpy.settings.preset` is :attr:`~scanpy.Preset.ScanpyV2Preview`,
        :class:`str`\ s are :meth:`anndata.acc.AdAcc.resolve`\ d (e.g. `'obsm.umap'`),
        otherwise they name an embedding found in `adata.obsm['X_[basis]']`
        (or `adata.obsm['[basis]']`).
        By default, the UMAP embedding.
    groupby
        Categorical observation/cell annotation for which densities
        are calculated per category.
        Given by a reference like `A.obs['phase']`, or by column name
        (:class:`str`\ s are :meth:`~anndata.acc.AdAcc.resolve`\ d with :attr:`~scanpy.Preset.ScanpyV2Preview`).
    components
        The embedding dimensions over which the density should be calculated.
        This is limited to two components.
    {out}
    basis
        The deprecated name of `use`.

    Returns
    -------
    Returns the densities if `out=True`, else `None`.
    Sets the following fields for `key` as `[basis]_density_[groupby]`, where `[basis]` is one of `umap`, `diffmap`, `pca`, `tsne`, or `draw_graph_fa` and `[groupby]` denotes the parameter input:

    `adata.obs[key]` : :class:`numpy.ndarray` (dtype `float`)
        Embedding density values for each cell.
    `adata.uns['[key]_params']` : :class:`dict`
        A dict with the values for the parameters `covariate` (for the `groupby` parameter) and `components`.

    Examples
    --------

    ..  exec-jupyter::

        import scanpy as sc
        adata = sc.datasets.pbmc68k_reduced()
        sc.tl.umap(adata)
        sc.tl.embedding_density(adata, 'umap', groupby='phase')
        sc.pl.embedding_density(
            adata, basis='umap', key='umap_density_phase', group='G1'
        )

    ..  exec-jupyter::

        sc.pl.embedding_density(
            adata, basis='umap', key='umap_density_phase', group='S'
        )

    .. currentmodule:: scanpy

    See Also
    --------
    pl.embedding_density

    """
    if basis is not None:
        if not isinstance(use, Default):
            msg = "Pass either `use` or `basis`, not both."
            raise TypeError(msg)
        use = basis
    basis = use
    groupby = _obs_col(groupby, param="groupby")
    if isinstance(basis, Default):
        basis = "umap"  # found in `X_umap` or `umap`
    elif isinstance(basis, str):
        basis = _resolve_ref(basis, vec=False)
    if isinstance(basis, MultiAcc):
        basis_key, basis = basis.k, basis.k.removeprefix("X_")
    else:
        basis = basis.lower()
        if basis == "fa":
            basis = "draw_graph_fa"
        basis_key = _get_basis_key(adata, basis)
    default = f"{basis}_density" + ("" if groupby is None else f"_{groupby}")
    key = _out_key(out, key_added, kind="obs", default=default)
    # to ensure that newly created covariates are categorical
    # to test for category numbers
    if key is not True:
        sanitize_anndata(adata)

    logg.info(f"computing density on {basis!r}")

    # Test user inputs
    if basis_key is None:
        msg = (
            "Cannot find the embedded representation "
            f"`adata.obsm[{basis!r} | 'X_{basis}']`. Compute the embedding first."
        )
        raise ValueError(msg)

    if components is None:
        components = "1,2"
    if isinstance(components, str):
        components = components.split(",")
    components = np.array(components).astype(int) - 1

    if len(components) != 2:
        msg = "Please specify exactly 2 components, or `None`."
        raise ValueError(msg)

    if basis == "diffmap":
        components += 1

    if groupby is not None:
        if groupby not in adata.obs:
            msg = f"Could not find {groupby!r} `.obs` column."
            raise ValueError(msg)

        groups = adata.obs[groupby]
        if key is True and pd.api.types.is_string_dtype(groups):
            groups = groups.astype("category")  # like `sanitize_anndata`, but local
        if groups.dtype.name != "category":
            msg = f"{groupby!r} column does not contain categorical data"
            raise ValueError(msg)

    # Calculate the densities over each category in the groupby column
    if groupby is not None:
        density_values = np.zeros(adata.n_obs)

        for cat in groups.cat.categories:
            cat_mask = (groups == cat).to_numpy()
            embed_x = adata.obsm[basis_key][cat_mask, components[0]]
            embed_y = adata.obsm[basis_key][cat_mask, components[1]]

            dens_embed = _calc_density(embed_x, embed_y)
            density_values[cat_mask] = dens_embed
    else:  # if groupby is None
        # Calculate the density over the whole embedding without subsetting
        embed_x = adata.obsm[basis_key][:, components[0]]
        embed_y = adata.obsm[basis_key][:, components[1]]

        density_values = _calc_density(embed_x, embed_y)
    if key is True:
        return density_values
    adata.obs[key] = density_values

    # Reduce diffmap components for labeling
    # Note: plot_scatter takes care of correcting diffmap components
    #       for plotting automatically
    if basis != "diffmap":
        components += 1

    adata.uns[f"{key}_params"] = dict(covariate=groupby, components=components.tolist())

    logg.hint(
        f"added\n"
        f"    '{key}', densities (adata.obs)\n"
        f"    '{key}_params', parameter (adata.uns)"
    )
    return None
