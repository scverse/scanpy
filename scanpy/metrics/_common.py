from __future__ import annotations

from functools import singledispatch
import warnings

import numpy as np
import pandas as pd
from scipy import sparse

from .._compat import DaskArray


@singledispatch
def _resolve_vals(val):
    return np.asarray(val)


@_resolve_vals.register(np.ndarray)
@_resolve_vals.register(sparse.csr_matrix)
@_resolve_vals.register(DaskArray)
def _(val):
    return val


@_resolve_vals.register(sparse.spmatrix)
def _(val):
    return sparse.csr_matrix(val)


@_resolve_vals.register(pd.DataFrame)
@_resolve_vals.register(pd.Series)
def _(val):
    return val.to_numpy()


def _check_vals(vals):
    """\
    Checks that values wont cause issues in computation.

    Returns new set of vals, and indexer to put values back into result.

    For details on why this is neccesary, see:
    https://github.com/scverse/scanpy/issues/1806
    """
    from scanpy._utils import is_constant

    full_result = np.empty(vals.shape[0], dtype=np.float64)
    full_result.fill(np.nan)
    idxer = ~is_constant(vals, axis=1)
    if idxer.all():
        idxer = slice(None)
    else:
        warnings.warn(
            UserWarning(
                f"{len(idxer) - idxer.sum()} variables were constant, will return nan for these.",
            )
        )
    return vals[idxer], idxer, full_result
