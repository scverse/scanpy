"""
Metrics which don't quite deserve their own file.
"""
from typing import Sequence, Union
import pandas as pd
import numpy as np


def confusion_matrix(
    orig: "Union[pd.Series, np.ndarray, Sequence]",
    new: "Union[pd.Series, np.ndarray, Sequence]",
    data: pd.DataFrame = None,
    *,
    normalize: bool = True,
) -> pd.DataFrame:
    """Given an original and new set of labels, create a labelled confusion matrix.

    Params
    ------
    orig
        Original labels
    new
        New labels
    normalize
        Should the confusion matrix be normalized?

    Usage
    -----
    >>> import seaborn as sns
    >>> import scanpy as sc
    >>> confmtx = sc.metrics.confusion_matrix(labels["orig"], labels["new"])
    >>> sns.heatmap(confmtx)
    """
    from sklearn.metrics import confusion_matrix as _confusion_matrix

    if data is not None:
        orig = data[orig]
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
    if orig.name is None:
        orig_name = "Original labels"
    else:
        orig_name = orig.name
    if new.name is None:
        new_name = "New labels"
    else:
        new_name = new.name
    df = pd.DataFrame(
        mtx,
        index=pd.Index(unique_labels, name=orig_name),
        columns=pd.Index(unique_labels, name=new_name),
    )

    # Filter
    df = df.loc[df.index.isin(pd.unique(orig)), df.columns.isin(pd.unique(new))]

    return df
