from typing import Optional

import numpy as np
from .. import logging as logg
from ._pca import pca
from .. import settings
from .._utils import _choose_graph

doc_use_rep = """\
use_rep
    Use the indicated representation. `'X'` or any key for `.obsm` is valid.
    If `None`, the representation is chosen automatically:
    For `.n_vars` < :attr:`scanpy.settings.N_PCS` (default: 50), `.X` is used, otherwise 'X_pca' is used.
    If 'X_pca' is not present, it’s computed with default parameters or `n_pcs` if present.\
"""

doc_n_pcs = """\
n_pcs
    Use this many PCs. If `n_pcs==0` use `.X` if `use_rep is None`.\
"""


def _choose_representation(adata, use_rep=None, n_pcs=None, silent=False):
    verbosity = settings.verbosity
    if silent and settings.verbosity > 1:
        settings.verbosity = 1
    if use_rep is None and n_pcs == 0:  # backwards compat for specifying `.X`
        use_rep = 'X'
    if use_rep is None:
        if adata.n_vars > settings.N_PCS:
            if 'X_pca' in adata.obsm.keys():
                if n_pcs is not None and n_pcs > adata.obsm['X_pca'].shape[1]:
                    raise ValueError(
                        '`X_pca` does not have enough PCs. Rerun `sc.pp.pca` with adjusted `n_comps`.'
                    )
                X = adata.obsm['X_pca'][:, :n_pcs]
                logg.info(f'    using \'X_pca\' with n_pcs = {X.shape[1]}')
            else:
                logg.warning(
                    f'You’re trying to run this on {adata.n_vars} dimensions of `.X`, '
                    'if you really want this, set `use_rep=\'X\'`.\n         '
                    'Falling back to preprocessing with `sc.pp.pca` and default params.'
                )
                n_pcs_pca = n_pcs if n_pcs is not None else settings.N_PCS
                X = pca(adata.X, n_comps=n_pcs_pca)
                adata.obsm['X_pca'] = X
        else:
            logg.info('    using data matrix X directly')
            X = adata.X
    else:
        if use_rep in adata.obsm.keys() and n_pcs is not None:
            if n_pcs > adata.obsm[use_rep].shape[1]:
                raise ValueError(
                    f'{use_rep} does not have enough Dimensions. Provide a '
                    'Representation with equal or more dimensions than'
                    '`n_pcs` or lower `n_pcs` '
                )
            X = adata.obsm[use_rep][:, :n_pcs]
        elif use_rep in adata.obsm.keys() and n_pcs is None:
            X = adata.obsm[use_rep]
        elif use_rep == 'X':
            X = adata.X
        else:
            raise ValueError(
                'Did not find {} in `.obsm.keys()`. '
                'You need to compute it first.'.format(use_rep)
            )
    settings.verbosity = verbosity  # resetting verbosity
    return X


def preprocess_with_pca(adata, n_pcs: Optional[int] = None, random_state=0):
    """
    Parameters
    ----------
    n_pcs
        If `n_pcs=0`, do not preprocess with PCA.
        If `None` and there is a PCA version of the data, use this.
        If an integer, compute the PCA.
    """
    if n_pcs == 0:
        logg.info('    using data matrix X directly (no PCA)')
        return adata.X
    elif n_pcs is None and 'X_pca' in adata.obsm_keys():
        logg.info(f'    using \'X_pca\' with n_pcs = {adata.obsm["X_pca"].shape[1]}')
        return adata.obsm['X_pca']
    elif 'X_pca' in adata.obsm_keys() and adata.obsm['X_pca'].shape[1] >= n_pcs:
        logg.info(f'    using \'X_pca\' with n_pcs = {n_pcs}')
        return adata.obsm['X_pca'][:, :n_pcs]
    else:
        n_pcs = settings.N_PCS if n_pcs is None else n_pcs
        if adata.X.shape[1] > n_pcs:
            logg.info(f'    computing \'X_pca\' with n_pcs = {n_pcs}')
            logg.hint('avoid this by setting n_pcs = 0')
            X = pca(adata.X, n_comps=n_pcs, random_state=random_state)
            adata.obsm['X_pca'] = X
            return X
        else:
            logg.info('    using data matrix X directly (no PCA)')
            return adata.X


def get_init_pos_from_paga(
    adata, adjacency=None, random_state=0, neighbors_key=None, obsp=None
):
    np.random.seed(random_state)
    if adjacency is None:
        adjacency = _choose_graph(adata, obsp, neighbors_key)
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
                init_pos[subset] = group_pos - 0.5 * dist + noise
            else:
                init_pos[subset] = group_pos
    else:
        raise ValueError(
            'Plot PAGA first, so that adata.uns[\'paga\']' 'with key \'pos\'.'
        )
    return init_pos
