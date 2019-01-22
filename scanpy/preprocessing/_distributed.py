import numpy as np

# install dask if available
try:
    import dask.array as da
except ImportError:
    da = None

def materialize_as_ndarray(a):
    """Convert distributed arrays to ndarrays."""
    if type(a) in (list, tuple):
        if da is not None and any(isinstance(arr, da.Array) for arr in a):
            return da.compute(*a, sync=True)
        return tuple(np.asarray(arr) for arr in a)
    return np.asarray(a)
