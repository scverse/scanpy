"""
Metrics which don't quite deserve their own file.
"""
from typing import Optional, Sequence, Union

import pandas as pd
from pandas.api.types import is_categorical_dtype
from natsort import natsorted
import numpy as np


def confusion_matrix(
    orig: Union[pd.Series, np.ndarray, Sequence],
    new: Union[pd.Series, np.ndarray, Sequence],
    data: Optional[pd.DataFrame] = None,
    *,
    normalize: bool = True,
) -> pd.DataFrame:
    """\
    Given an original and new set of labels, create a labelled confusion matrix.

    Parameters `orig` and `new` can either be entries in data or categorical arrays
    of the same size.

    Params
    ------
    orig
        Original labels.
    new
        New labels.
    data
        Optional dataframe to fill entries from.
    normalize
        Should the confusion matrix be normalized?


    Examples
    --------

    .. plot::

        import scanpy as sc; import seaborn as sns
        pbmc = sc.datasets.pbmc68k_reduced()
        cmtx = sc.metrics.confusion_matrix("bulk_labels", "louvain", pbmc.obs)
        sns.heatmap(cmtx)

    """
    from sklearn.metrics import confusion_matrix as _confusion_matrix

    if data is not None:
        if isinstance(orig, str):
            orig = data[orig]
        if isinstance(new, str):
            new = data[new]

    # Coercing so I don't have to deal with it later
    orig, new = pd.Series(orig), pd.Series(new)
    assert len(orig) == len(new)

    unique_labels = pd.unique(np.concatenate((orig.values, new.values)))

    # Compute
    mtx = _confusion_matrix(orig, new, labels=unique_labels)
    if normalize:
        sums = mtx.sum(axis=1)[:, np.newaxis]
        mtx = np.divide(mtx, sums, where=sums != 0)

    # Label
    orig_name = "Original labels" if orig.name is None else orig.name
    new_name = "New Labels" if new.name is None else new.name
    df = pd.DataFrame(
        mtx,
        index=pd.Index(unique_labels, name=orig_name),
        columns=pd.Index(unique_labels, name=new_name),
    )

    # Filter
    if is_categorical_dtype(orig):
        orig_idx = pd.Series(orig).cat.categories
    else:
        orig_idx = natsorted(pd.unique(orig))
    if is_categorical_dtype(new):
        new_idx = pd.Series(new).cat.categories
    else:
        new_idx = natsorted(pd.unique(new))
    df = df.loc[np.array(orig_idx), np.array(new_idx)]

    return df
