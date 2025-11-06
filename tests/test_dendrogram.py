from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

import scanpy as sc
from testing.scanpy._helpers.data import pbmc68k_reduced

n_neighbors = 5
key = "test"


@pytest.mark.parametrize("groupby", ["bulk_labels", ["bulk_labels", "phase"]])
@pytest.mark.parametrize("key_added", [None, "custom_key"])
def test_dendrogram_key_added(groupby, key_added):
    adata = pbmc68k_reduced()
    sc.tl.dendrogram(adata, groupby=groupby, key_added=key_added, use_rep="X_pca")
    if isinstance(groupby, list):
        dendrogram_key = f"dendrogram_{'_'.join(groupby)}"
    else:
        dendrogram_key = f"dendrogram_{groupby}"

    if key_added is None:
        key_added = dendrogram_key
    assert key_added in adata.uns


REP_PCA_0 = [
    *(1.50808525e00, -1.67258829e-01, -7.12063432e-01, -2.07935140e-01),
    *(-3.55730444e-01, -2.24421427e-01, -1.46907698e-02, -7.01090470e-02),
    *(-1.31467551e-01, -3.75757217e-02, -1.07698059e-02, -4.37555499e-02),
    *(1.06897885e-02, 1.10454357e-03, -5.37674241e-02, -4.94170748e-03),
    *(1.11988001e-02, -4.48330259e-03, -2.56892946e-02, -3.50749046e-02),
    *(-3.15931924e-02, 2.84416862e-02, -3.70664597e-02, -2.38820408e-02),
    *(-4.57040370e-02, 2.99325008e-02, 9.56365839e-03, -4.28026691e-02),
    *(5.36734704e-03, -3.08445804e-02, -1.16719725e-02, -2.35078149e-02),
    *(2.87542702e-04, -1.70532353e-02, -1.79676879e-02, -3.09410989e-02),
    *(-1.09178647e-02, -1.60753895e-02, -1.04769412e-02, -1.36501975e-02),
    *(-6.83976896e-03, -1.17562497e-02, -4.65345643e-02, 1.91588048e-02),
    *(-1.38043752e-02, 4.75460896e-03, -1.41307563e-02, -1.03387292e-02),
    *(-1.68043356e-02, 1.33516011e-03),
]
REP_PCA_1_6 = [
    *(-2.70745814e-01, -3.45929652e-01, 6.27844110e-02, -8.34012777e-02),
    *(-1.08290315e-01, -1.38125733e-01, -2.57148240e-02, -2.73127705e-02),
    *(-1.45030200e-01, -6.88858554e-02, -4.28490154e-02, -1.88931823e-02),
    *(-2.56232135e-02, -7.66322482e-03, -5.49384989e-02, -1.43514248e-02),
    *(2.42769364e-02, -3.01547404e-02, -3.37253511e-02, -3.81337740e-02),
    *(-3.42049589e-03, -4.34436463e-03, -4.15385924e-02, -2.66448390e-02),
    *(-2.74285320e-02, 1.47806173e-02, 1.19129466e-02, -6.70884028e-02),
    *(2.58150720e-03, -1.64280720e-02, -1.07431635e-02, -3.04328315e-02),
    *(-3.82748269e-03, -2.95090005e-02, -3.10521629e-02, -3.43420058e-02),
    *(-4.49432433e-03, -2.15906072e-02, -1.23507539e-02, -2.88041346e-02),
    *(-7.31994957e-03, -7.28111062e-03, -7.61008039e-02, 2.40524579e-02),
    *(-1.20806806e-02, 5.05997473e-03, -2.53410172e-02, -1.83318909e-02),
    *(-1.81263424e-02, -3.35110351e-03),
]
REP_PCA = np.array([REP_PCA_0, *([REP_PCA_1_6] * 6)], dtype=np.float32)


def test_dendrogram_cor():
    rep = sc.AnnData(
        sparse.csr_matrix(  # noqa: TID251
            (
                np.array([1.2762934659055623, 1.6916760106710726, 1.6916760106710726]),
                np.array([12, 5, 44]),
                np.array([0, 0, 0, 0, 0, 1, 2, 3]),
            ),
            shape=(7, 51),
        ),
        dict(leiden=pd.Categorical(["372", "366", "357", "357", "357", "357", "357"])),
        obsm=dict(X_pca=REP_PCA),
    )
    sc.tl.dendrogram(rep, groupby="leiden")
