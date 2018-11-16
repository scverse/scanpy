# install dask if available
try:
    import dask.array as da
except ImportError:
    da = None

# install zappy (which wraps numpy), or fall back to plain numpy
try:
    import zappy.base as np
except ImportError:
    import numpy as np


def materialize_as_ndarray(a):
    """Convert distributed arrays to ndarrays."""
    if type(a) in (list, tuple):
        if da is not None and any(isinstance(arr, da.Array) for arr in a):
            return da.compute(*a, sync=True)
        elif hasattr(np, 'asarray'): # zappy case
            return tuple(np.asarray(arr) for arr in a)
    else:
        if da is not None and isinstance(a, da.Array):
            return a.compute()
        elif hasattr(np, 'asarray'): # zappy case
            return np.asarray(a)
    return a
