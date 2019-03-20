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

    return(z)


def density(
        adata: AnnData,
        embedding: str,
        key: str = None):
    """Calculate the density of cells in an embedding (per condition)

    Gaussian kernel density estimation is used to calculate the density of
    cells in an embedded space. This can be performed per category over a
    categorical cell annotation. The cell density can be plotted using the 
    `sc.pl.density()` function.

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

    Returns
    -------
    Updates `adata.obs` with an additional field `[embedding]_density`, where
    `[embedding]` is one of `umap`, `dm`, `pca`, `tsne`, or `draw_graph_fa`.
    Updates `adata.uns` with an additional field `[embedding]_density_param`.
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

    # Calculate the densities over each category in the groupby column
    if groupby is not None:
        categories = adata.obs[groupby].cat.categories

        adata.obs[embedding+'_density'] = [0 for i in range(adata.n_obs)]
        
        for cat in categories:
            cat_mask = adata.obs[groupby] == cat
            embed_x = adata.obsm['X_'+embedding][cat_mask, 0]
            embed_y = adata.obsm['X_'+embedding][cat_mask, 1]

            dens_embed = _calc_density(embed_x, embed_y)
            adata.obs[embedding+'_density'][cat_mask] = dens_embed

    # Calculate the density over the whole embedding without subsetting
    else: #if groupby is None
        embed_x = adata.obsm['X_'+embedding][:, 0]
        embed_y = adata.obsm['X_'+embedding][:, 1]

        adata.obs[embedding+'_density'] = _calc_density(embed_x, embed_y)

    adata.uns[embedding+'_density_param'] = groupby

    logg.hint('added\n'
              '    \'{}_density\', densities (adata.obs)\n'
              '    \'{}_density_param\', parameter (adata.uns)'.format(embedding, embedding))

    return None
    


    
