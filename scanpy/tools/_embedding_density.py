"""Calculate density of cells in embeddings
"""

import numpy as np
from scipy.stats import gaussian_kde
from anndata import AnnData
from typing import Union, Optional, Sequence

from .. import logging as logg
from ..utils import sanitize_anndata

def _calc_density(
    x: np.ndarray,
    y: np.ndarray,
):
    """
    Function to calculate the density of cells in an embedding.
    """    

    # Calculate the point density
    xy = np.vstack([x,y])
    z = gaussian_kde(xy)(xy)

    min_z = np.min(z)
    max_z = np.max(z)

    # Scale between 0 and 1
    scaled_z = (z-min_z)/(max_z-min_z)

    return(scaled_z)


def embedding_density(
    adata: AnnData,
    basis: str,
    *,
    groupby: Optional[str] = None,
    key_added: Optional[str] = None,
    components: Union[str, Sequence[str]] = None
):
    """Calculate the density of cells in an embedding (per condition)

    Gaussian kernel density estimation is used to calculate the density of
    cells in an embedded space. This can be performed per category over a
    categorical cell annotation. The cell density can be plotted using the 
    `sc.pl.embedding_density()` function.

    Note that density values are scaled to be between 0 and 1. Thus, the
    density value at each cell is only comparable to other densities in 
    the same condition category.

    This function was written by Sophie Tritschler and implemented into
    Scanpy by Malte Luecken.
    
    Parameters
    ----------
    adata
        The annotated data matrix.
    basis
        The embedding over which the density will be calculated. This embedded
        representation should be found in `adata.obsm['X_[basis]']``.
    groupby
        Keys for categorical observation/cell annotation for which densities
        are calculated per category. Columns with up to ten categories are
        accepted.
    key_added
        Name of the `.obs` covariate that will be added with the density
        estimates.
    components
        The embedding dimensions over which the density should be calculated.
        This is limited to two components.

    Returns
    -------
    Updates `adata.obs` with an additional field specified by the `key_added`
    parameter. This parameter defaults to `[basis]_density_[groupby]`, where
    where `[basis]` is one of `umap`, `diffmap`, `pca`, `tsne`, or `draw_graph_fa`
    and `[groupby]` denotes the parameter input.
    Updates `adata.uns` with an additional field `[key_added]_params`.

    Examples
    --------
    >>> adata = sc.datasets.pbmc68k_reduced()
    >>> sc.tl.umap(adata)
    >>> sc.tl.embedding_density(adata, basis='umap', groupby='phase')
    >>> sc.pl.embedding_density(adata, basis='umap', key='umap_density_phase', 
    ...                         group='G1')
    >>> sc.pl.embedding_density(adata, basis='umap', key='umap_density_phase', 
    ...                         group='S')
    """
    sanitize_anndata(adata) # to ensure that newly created covariates are categorical to test for categoy numbers

    logg.info('computing density on \'{}\''.format(basis), r=True)

    # Test user inputs
    basis = basis.lower()
    
    if basis == 'fa':
        basis = 'draw_graph_fa'

    if 'X_'+basis not in adata.obsm_keys():
        raise ValueError('Cannot find the embedded representation `adata.obsm[X_{!r}]`. '
                         'Compute the embedding first.'.format(basis))

    if components is None: components = '1,2'
    if isinstance(components, str): components = components.split(',')
    components = np.array(components).astype(int) - 1

    if len(components) != 2:
        raise ValueError('Please specify exactly 2 components, or `None`.')

    if basis == 'diffmap': components += 1

    if groupby is not None:
        if groupby not in adata.obs:
            raise ValueError('Could not find {!r} `.obs` column.'.format(groupby))

        if adata.obs[groupby].dtype.name != 'category':
            raise ValueError('{!r} column does not contain Categorical data'.format(groupby))
    
        if len(adata.obs[groupby].cat.categories) > 10:
            raise ValueError('More than 10 categories in {!r} column.'.format(groupby))

    # Define new covariate name
    if key_added is not None:
        density_covariate = key_added
    elif groupby is not None:
        density_covariate = basis+'_density_'+groupby
    else:
        density_covariate = basis+'_density'

    # Calculate the densities over each category in the groupby column
    if groupby is not None:
        categories = adata.obs[groupby].cat.categories

        density_values = np.zeros(adata.n_obs)
        
        for cat in categories:
            cat_mask = adata.obs[groupby] == cat
            embed_x = adata.obsm['X_'+basis][cat_mask, components[0]]
            embed_y = adata.obsm['X_'+basis][cat_mask, components[1]]

            dens_embed = _calc_density(embed_x, embed_y)
            density_values[cat_mask] = dens_embed

        adata.obs[density_covariate] = density_values
        
    # Calculate the density over the whole embedding without subsetting
    else: #if groupby is None
        embed_x = adata.obsm['X_'+basis][:, components[0]]
        embed_y = adata.obsm['X_'+basis][:, components[1]]

        adata.obs[density_covariate] = _calc_density(embed_x, embed_y)


    # Reduce diffmap components for labeling
    # Note: plot_scatter takes care of correcting diffmap components for plotting automatically
    if basis != 'diffmap': components += 1

    adata.uns[density_covariate+'_params'] = {'covariate':groupby, 'components':components.tolist()}
    

    logg.hint('added\n'
              '    \'{}\', densities (adata.obs)\n'
              '    \'{}_params\', parameter (adata.uns)'.format(density_covariate, density_covariate))

    return None
