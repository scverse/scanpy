from typing import Optional

import numpy as np
from scipy import sparse
from sklearn.utils import check_random_state

from anndata import AnnData


def _choose_obs_rep(adata, *, use_raw=False, layer=None, obsm=None, obsp=None):
    """
    Choose array aligned with obs annotation.
    """
    is_layer = layer is not None
    is_raw = use_raw is not False
    is_obsm = obsm is not None
    is_obsp = obsp is not None
    choices_made = sum((is_layer, is_raw, is_obsm, is_obsp))
    assert choices_made <= 1
    if choices_made == 0:
        return adata.X
    elif is_layer:
        return adata.layers[layer]
    elif use_raw:
        return adata.raw.X
    elif is_obsm:
        return adata.obsm[obsm]
    elif is_obsp:
        return adata.obsp[obsp]
    else:
        raise RuntimeError("You broke it. But how? Please report this.")


def ica(
    adata: AnnData,
    n_components: int,
    *,
    random_state: Optional[int] = None,
    whiten: bool = True,
    mean_center: bool = True,
    inplace: bool = True,
    use_raw: bool = False,
    layer: Optional[str] = None,
    obsm: Optional[str] = None,
    highly_variable_genes: bool = False,
    **kwargs
) -> "Optional[FastICA]":
    """
    Independent components analysis.

    Computes observation and variable loadings of an ICA decomposition. This is
    based on (and uses parts of) the scikit-learn implementation (see
    :class:`~sklearn.decomposition.FastICA` but with a more efficient whitening
    implementation.

    Params
    ------
    # TODO

    Usage
    -----

    >>> pbmc = sc.datasets.pbmc64k_reduced()
    >>> sc.tl.ica(pbmc, 10, use_raw=True)
    # Using ICA coords for nearest neigbors and clustering:
    >>> sc.tl.neighbors(pbmc, use_rep="X_ica")  
    >>> sc.tl.leiden(pbmc)
    # Plotting ICA coords
    >>> sc.pl.embedding(pbmc, basis="ica", color="leiden")
    """
    from sklearn.decomposition import FastICA
    random_state = check_random_state(random_state)

    X = _choose_obs_rep(adata, use_raw=use_raw, layer=layer, obsm=obsm)
    if highly_variable_genes:
        X = X[:, adata.var["highly_variable"].values].copy()
    X = X.T

    ica_transformer = FastICA(n_components=n_components, whiten=False, **kwargs)

    if whiten:
        # Code heavily based on sklearn FastICA implementation
        # Swapped out np.lingalg.svd with scipy.sparse.linalg.svds
        # TODO: Make sure this is actually whitening (https://github.com/scikit-learn/scikit-learn/pull/13096)
        # TODO: Make sure what I do is equivalent to what sklearn does
        # Mean center:
        if mean_center:
            if sparse.issparse(X):
                X = X.toarray()
            else:
                X = X.copy()
            X_mean = X.mean(axis=-1)
            X -= X_mean[:, np.newaxis]
            del X_mean
        # Using arpack for whitening, could add other options as well
        # Maybe use sklearn transformer for this?
        v0 = random_state.uniform(-1, 1, size=min(X.shape))
        u, d, _ = sparse.linalg.svds(X, k=n_components, v0=v0)
        # svds doesn't abide by scipy.linalg.svd/randomized_svd
        # conventions, so reverse its outputs.
        d = d[::-1]
        K = (u / d).T[:n_components]
        X1 = K @ X  # np.dot takes way more memory if X is sparse
        X1 *= np.sqrt(X.shape[1])

        del u, _, v0
    else:
        raise NotImplementedError()

    ica_transformer.fit(X1.T)

    # Handle transformation from whitening
    if whiten:
        W = ica_transformer.components_
        S = (W @ K @ X).T
        components = W @ K

    else:
        components = ica_transformer.components_
        S = (components @ X).T

    if inplace:
        adata.obsm["X_ica"] = S
        if highly_variable_genes:
            ICs = np.zeros((adata.n_vars, n_components))
            ICs[adata.var["highly_variable"].values, :] = components.T
        else:
            ICs = components.T
        adata.varm["ICs"] = ICs
    else:
        return ica_transformer