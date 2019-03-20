"""Calculate density of cells in embeddings
"""

import numpy as np
from scipy.stats import gaussian_kde
from anndata import AnnData

from .. import logging as logg


def _calc_density(
        x: np.ndarray,
        y: np.ndarray):
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


def density(
        adata: AnnData,
        embedding: str,
        groupby: str = None,
        key_added: str = None):
    """Calculate the density of cells in an embedding (per condition)

    Gaussian kernel density estimation is used to calculate the density of
    cells in an embedded space. This can be performed per category over a
    categorical cell annotation. The cell density can be plotted using the 
    `sc.pl.density()` function.

    Note that density values are scaled to be between 0 and 1, so that
    the density at each cell is only comparable to other densities in the
    same condition category.

    This function was written by Sophie Tritschler and implemented into
    Scanpy by Malte Luecken.
    
    Parameters
    ----------
    adata : :class:`~anndata.AnnData`
        The annotated data matrix.
    embedding : `str`
        The embedding over which the density will be calculated. This must
        be one of:
        'umap' : UMAP
        'dm' : Diffusion map
        'pca' : PCA
        'tsne' : t-SNE
        'fa' : Force-directed graph layout by Force Atlas 2
    groupby : `str`, optional (default: `None`)
        Keys for categorical observation/cell annotation for which densities
        are calculated per category. Columns with up to ten categories are
        accepted.
    key_added : `str`, optional (default: `None`)
        Name of the `.obs` covariate that will be added with the density
        estimates.

    Returns
    -------
    Updates `adata.obs` with an additional field specified by the `key_added`
    parameter. This parameter defaults to `[embedding]_density_[groupby]`,
    where `[embedding]` is one of `umap`, `dm`, `pca`, `tsne`, or `draw_graph_fa`
    and `[groupby]` denotes the parameter input.
    Updates `adata.uns` with an additional field `[key_added]_param`.
    """
    logg.info('computing density on \'{}\''.format(embedding), r=True)

    # Test if inputs are okay
    embedding = embedding.lower()
    
    if embedding == 'fa':
        embedding = 'draw_graph_fa'

    allowed_embeddings = ['umap', 'dm', 'pca', 'tsne', 'draw_graph_fa']

    if embedding not in allowed_embeddings:
        raise ValueError('{!r} is not a valid embedding.'.format(embedding))

    if 'X_'+embedding not in adata.obsm:
        raise ValueError('Cannot find the embedded representation. Compute the embedding first.')

    if groupby is not None:
        if groupby not in adata.obs:
            raise ValueError('Could not find {!r} `.obs` column.'.format(groupby))

        if adata.obs[groupby].dtype.name == 'category':
            raise ValueError('{!r} column does not contain Categorical data'.format(groupby))
    
        if len(adata.obs[groupby].cat.categories) > 10:
            raise ValueError('More than 10 categories in {!r} column.'.format(groupby))

    # Define new covariate name
    if key_added is not None:
        density_covariate = key_added
    elif groupby is not None:
        density_covariate = embedding+'_density_'+groupby
    else:
        density_covariate = embedding+'_density'

    # Calculate the densities over each category in the groupby column
    if groupby is not None:
        categories = adata.obs[groupby].cat.categories

        adata.obs[density_covariate] = [0 for i in range(adata.n_obs)]
        
        for cat in categories:
            cat_mask = adata.obs[groupby] == cat
            embed_x = adata.obsm['X_'+embedding][cat_mask, 0]
            embed_y = adata.obsm['X_'+embedding][cat_mask, 1]

            dens_embed = _calc_density(embed_x, embed_y)
            adata.obs[density_covariate][cat_mask] = dens_embed

    # Calculate the density over the whole embedding without subsetting
    else: #if groupby is None
        embed_x = adata.obsm['X_'+embedding][:, 0]
        embed_y = adata.obsm['X_'+embedding][:, 1]

        adata.obs[density_covariate] = _calc_density(embed_x, embed_y)

    adata.uns[density_covariate+'_param'] = groupby

    logg.hint('added\n'
              '    \'{}\', densities (adata.obs)\n'
              '    \'{}_param\', parameter (adata.uns)'.format(density_covariate, density_covariate))

    return None
    


    
