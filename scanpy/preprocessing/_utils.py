import numpy as np
from scipy.sparse import issparse


def _get_mean_var(X, weights=None):
    # - using sklearn.StandardScaler throws an error related to
    #   int to long trafo for very large matrices
    # - using X.multiply is slower
    if True:
        if not issparse(X):
            mean = np.average(X, axis=0, weights=weights)
            mean_sq = np.average(X.power(2), axis=0, weights=weights)
        elif weights is not None:
            mean = np.average(X.toarray(), axis=0, weights=weights)
            mean_sq = np.average(X.power(2).toarray(), axis=0, weights=weights)
        else:
            mean = X.mean(axis=0).A1
            mean_sq = X.power(2).mean(axis=0).A1
        # enforece R convention (unbiased estimator) for variance
        var = (mean_sq - mean**2) * (X.shape[0]/(X.shape[0]-1))
    else:
        from sklearn.preprocessing import StandardScaler
        scaler = StandardScaler(with_mean=False).partial_fit(X)
        mean = scaler.mean_
        # enforce R convention (unbiased estimator)
        var = scaler.var_ * (X.shape[0]/(X.shape[0]-1))
    return mean, var
