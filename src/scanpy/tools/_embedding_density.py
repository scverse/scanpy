"""Calculate density of cells in embeddings."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from .. import logging as logg
from .._compat import old_positionals
from .._utils import sanitize_anndata

if TYPE_CHECKING:
    from collections.abc import Sequence

    from anndata import AnnData


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


@old_positionals("groupby", "key_added", "components")
def embedding_density(  # noqa: PLR0912
    adata: AnnData,
    basis: str = "umap",
    *,
    groupby: str | None = None,
    key_added: str | None = None,
    components: str | Sequence[str] | None = None,
) -> None:
    """Calculate the density of cells in an embedding (per condition).

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

    Parameters
    ----------
    adata
        The annotated data matrix.
    basis
        The embedding over which the density will be calculated. This embedded
        representation is found in `adata.obsm['X_[basis]']``.
    groupby
        Key for categorical observation/cell annotation for which densities
        are calculated per category.
    key_added
        Name of the `.obs` covariate that will be added with the density
        estimates.
    components
        The embedding dimensions over which the density should be calculated.
        This is limited to two components.

    Returns
    -------
    Sets the following fields (`key_added` defaults to `[basis]_density_[groupby]`, where `[basis]` is one of `umap`, `diffmap`, `pca`, `tsne`, or `draw_graph_fa` and `[groupby]` denotes the parameter input):

    `adata.obs[key_added]` : :class:`numpy.ndarray` (dtype `float`)
        Embedding density values for each cell.
    `adata.uns['[key_added]_params']` : :class:`dict`
        A dict with the values for the parameters `covariate` (for the `groupby` parameter) and `components`.

    Examples
    --------

    .. plot::
        :context: close-figs

        import scanpy as sc
        adata = sc.datasets.pbmc68k_reduced()
        sc.tl.umap(adata)
        sc.tl.embedding_density(adata, basis='umap', groupby='phase')
        sc.pl.embedding_density(
            adata, basis='umap', key='umap_density_phase', group='G1'
        )

    .. plot::
        :context: close-figs

        sc.pl.embedding_density(
            adata, basis='umap', key='umap_density_phase', group='S'
        )

    .. currentmodule:: scanpy

    See Also
    --------
    pl.embedding_density

    """
    # to ensure that newly created covariates are categorical
    # to test for category numbers
    sanitize_anndata(adata)

    logg.info(f"computing density on {basis!r}")

    # Test user inputs
    basis = basis.lower()

    if basis == "fa":
        basis = "draw_graph_fa"

    if f"X_{basis}" not in adata.obsm_keys():
        msg = (
            "Cannot find the embedded representation "
            f"`adata.obsm['X_{basis}']`. Compute the embedding first."
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

        if adata.obs[groupby].dtype.name != "category":
            msg = f"{groupby!r} column does not contain categorical data"
            raise ValueError(msg)

    # Define new covariate name
    if key_added is not None:
        density_covariate = key_added
    elif groupby is not None:
        density_covariate = f"{basis}_density_{groupby}"
    else:
        density_covariate = f"{basis}_density"

    # Calculate the densities over each category in the groupby column
    if groupby is not None:
        categories = adata.obs[groupby].cat.categories

        density_values = np.zeros(adata.n_obs)

        for cat in categories:
            cat_mask = adata.obs[groupby] == cat
            embed_x = adata.obsm[f"X_{basis}"][cat_mask, components[0]]
            embed_y = adata.obsm[f"X_{basis}"][cat_mask, components[1]]

            dens_embed = _calc_density(embed_x, embed_y)
            density_values[cat_mask] = dens_embed

        adata.obs[density_covariate] = density_values
    else:  # if groupby is None
        # Calculate the density over the whole embedding without subsetting
        embed_x = adata.obsm[f"X_{basis}"][:, components[0]]
        embed_y = adata.obsm[f"X_{basis}"][:, components[1]]

        adata.obs[density_covariate] = _calc_density(embed_x, embed_y)

    # Reduce diffmap components for labeling
    # Note: plot_scatter takes care of correcting diffmap components
    #       for plotting automatically
    if basis != "diffmap":
        components += 1

    adata.uns[f"{density_covariate}_params"] = dict(
        covariate=groupby, components=components.tolist()
    )

    logg.hint(
        f"added\n"
        f"    '{density_covariate}', densities (adata.obs)\n"
        f"    '{density_covariate}_params', parameter (adata.uns)"
    )
