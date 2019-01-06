import numpy as np
from .. import logging as logg
from ._pca import pca
from ..preprocessing._simple import N_PCS

doc_use_rep = """\
use_rep : {`None`, 'X'} or any key for `.obsm`, optional (default: `None`)
    Use the indicated representation. If `None`, the representation is chosen
    automatically: for `.n_vars` < 50, `.X` is used, otherwise 'X_pca' is used.
    If 'X_pca' is not present, it's computed with default parameters.\
"""

doc_n_pcs = """\
n_pcs : `int` or `None`, optional (default: `None`)
    Use this many PCs. If `n_pcs==0` use `.X` if `use_rep is None`.\
"""


def choose_representation(adata, use_rep=None, n_pcs=None):
    if use_rep is None and n_pcs == 0:  # backwards compat for specifying `.X`
        use_rep = 'X'
    if use_rep is None:
        if adata.n_vars > N_PCS:
            if 'X_pca' in adata.obsm.keys():
                if n_pcs is not None and n_pcs > adata.obsm['X_pca'].shape[1]:
                    raise ValueError(
                        '`X_pca` does not have enough PCs. Rerun `sc.pp.pca` with adjusted `n_comps`.')
                X = adata.obsm['X_pca'][:, :n_pcs]
                logg.info('    using \'X_pca\' with n_pcs = {}'
                          .format(X.shape[1]))
                return X
            else:
                logg.warn(
                    'You\'re trying to run this on {} dimensions of `.X`, '
                    'if you really want this, set `use_rep=\'X\'`.\n         '
                    'Falling back to preprocessing with `sc.pp.pca` and default params.'
                    .format(adata.n_vars))
                X = pca(adata.X)
                adata.obsm['X_pca'] = X[:, :n_pcs]
                return X
        else:
            logg.info('    using data matrix X directly')
            return adata.X
    else:
        if use_rep in adata.obsm.keys():
            return adata.obsm[use_rep]
        elif use_rep == 'X':
            return adata.X
        else:
            raise ValueError(
                'Did not find {} in `.obsm.keys()`. '
                'You need to compute it first.'.format(use_rep))


def preprocess_with_pca(adata, n_pcs=None, random_state=0):
    """
    Parameters
    ----------
    n_pcs : `int` or `None`, optional (default: `None`)
        If `n_pcs=0`, do not preprocess with PCA.
        If `None` and there is a PCA version of the data, use this.
        If an integer, compute the PCA.
    """
    if n_pcs == 0:
        logg.info('    using data matrix X directly (no PCA)')
        return adata.X
    elif n_pcs is None and 'X_pca' in adata.obsm_keys():
        logg.info('    using \'X_pca\' with n_pcs = {}'
                  .format(adata.obsm['X_pca'].shape[1]))
        return adata.obsm['X_pca']
    elif ('X_pca' in adata.obsm_keys()
          and adata.obsm['X_pca'].shape[1] >= n_pcs):
        logg.info('    using \'X_pca\' with n_pcs = {}'
                  .format(n_pcs))
        return adata.obsm['X_pca'][:, :n_pcs]
    else:
        n_pcs = N_PCS if n_pcs is None else n_pcs
        if adata.X.shape[1] > n_pcs:
            logg.info('    computing \'X_pca\' with n_pcs = {}'.format(n_pcs))
            logg.hint('avoid this by setting n_pcs = 0')
            X = pca(adata.X, n_comps=n_pcs, random_state=random_state)
            adata.obsm['X_pca'] = X
            return X
        else:
            logg.info('    using data matrix X directly (no PCA)')
            return adata.X


def get_init_pos_from_paga(adata, adjacency=None, random_state=0):
    np.random.seed(random_state)
    if adjacency is None:
        adjacency = adata.uns['neighbors']['connectivities']
    if 'paga' in adata.uns and 'pos' in adata.uns['paga']:
        groups = adata.obs[adata.uns['paga']['groups']]
        pos = adata.uns['paga']['pos']
        connectivities_coarse = adata.uns['paga']['connectivities']
        init_pos = np.ones((adjacency.shape[0], 2))
        for i, group_pos in enumerate(pos):
            subset = (groups == groups.cat.categories[i]).values
            neighbors = connectivities_coarse[i].nonzero()
            if len(neighbors[1]) > 0:
                connectivities = connectivities_coarse[i][neighbors]
                nearest_neighbor = neighbors[1][np.argmax(connectivities)]
                noise = np.random.random((len(subset[subset]), 2))
                dist = pos[i] - pos[nearest_neighbor]
                noise = noise * dist
                init_pos[subset] = group_pos - 0.5*dist + noise
            else:
                init_pos[subset] = group_pos
    else:
        raise ValueError('Plot PAGA first, so that adata.uns[\'paga\']'
                         'with key \'pos\'.')
    return init_pos
