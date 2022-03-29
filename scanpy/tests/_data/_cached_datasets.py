from functools import wraps
import scanpy as sc


def cached_dataset(func):
    store = []

    @wraps(func)
    def wrapper():
        if len(store) < 1:
            store.append(func())
        return store[0].copy()

    return wrapper


pbmc3k = cached_dataset(sc.datasets.pbmc3k)
pbmc68k_reduced = cached_dataset(sc.datasets.pbmc68k_reduced)
pbmc3k_processed = cached_dataset(sc.datasets.pbmc3k_processed)
krumsiek11 = cached_dataset(sc.datasets.krumsiek11)
paul15 = cached_dataset(sc.datasets.paul15)
