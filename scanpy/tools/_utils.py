from .. import logging as logg
from .pca import pca
from ..preprocessing.simple import N_PCS
from textwrap import dedent

doc_use_rep = dedent("""\
    use_rep : \{`None`, 'X'\} or any key for `.obsm`, optional (default: `None`)
        Use the indicated representation. If `None`, the representation is
        chosen automatically: for `.n_vars` < 50, `.X` is used, otherwise and if
        'X_pca' is present, 'X_pca' is used.
""")


def choose_representation(adata, use_rep=None):
    if use_rep is None:
        if adata.n_vars > N_PCS:
            if 'X_pca' in adata.obsm.keys():
                logg.info('    using \'X_pca\' with n_pcs = {}'
                          .format(adata.obsm['X_pca'].shape[1]))
                return adata.obsm['X_pca']
            else:
                raise ValueError(
                    'You\`re trying to run the computation on `.n_vars` dimensions of `.X`.'
                    'If you really want this set `use_rep=\'X\'`. '
                    'Otherwise, for instance, perform preprocessing with PCA by calling `sc.pp.pca(adata)`')
        else:
            logg.info('    using data matrix X directly')
            return adata.X
    else:
        if use_rep in adata.obsm.keys():
            return adata.obsm[use_rep]
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
