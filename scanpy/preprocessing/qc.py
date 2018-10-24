import numpy as np
from scipy import sparse
import numba

def top_proportions(mtx, n):
    """Calculates cumulative proportions up to value n"""
    if sparse.issparse(mtx):
        if not sparse.isspmatrix_csr(mtx):
            mtx = sparse.csr_matrix(mtx)
        # Allowing numba to do more
        return top_proportions_sparse_csr(mtx.data, mtx.indptr, n)
    else:
        return top_proportions_dense(mtx, n)

@numba.jit
def top_proportions_dense(mtx, n):
    sums = mtx.sum(axis=1)
    partitioned = np.apply_along_axis(np.argpartition, 1, -mtx, n-1)
    partitioned = partitioned[:, :n]
    values = np.zeros_like(partitioned, dtype=np.float64)
    for i in range(partitioned.shape[0]):
        vec = mtx[i, partitioned[i, :]]  # Not a view
        vec[::-1].sort()  # Sorting on a reversed view (e.g. a descending sort)
        vec = np.cumsum(vec) / sums[i]
        values[i, :] = vec
    return values

@numba.jit(parallel=True)
def top_proportions_sparse_csr(data, indptr, n):
    values = np.zeros((indptr.size-1, n), dtype=np.float64)
    for i in numba.prange(indptr.size-1):
        start, end = indptr[i], indptr[i+1]
        vec = np.zeros(n, dtype=np.float64)
        if end - start <= n:
            vec[:end-start] = data[start:end]
            total = vec.sum()
        else:
            vec[:] = -(np.partition(-data[start:end], n-1)[:n])
            total = (data[start:end]).sum()  # Is this not just vec.sum()?
        vec[::-1].sort()
        values[i, :] = vec.cumsum() / total
    return values


def top_segment_proportions(mtx, ns):
    """Calculates total percentage of counts in top ns genes"""
    # Pretty much just does dispatch
    if sparse.issparse(mtx):
        if not sparse.isspmatrix_csr(mtx):
            mtx = sparse.csr_matrix(mtx)
        return top_segment_proportions_sparse_csr(mtx.data, mtx.indptr, ns)
    else:
        return top_segment_proportions_dense(mtx, ns)

def top_segment_proportions_dense(mtx, ns):
    # Currently ns is considered to be 1 indexed
    ns = np.sort(ns)
    sums = mtx.sum(axis=1)
    partitioned = np.apply_along_axis(
        np.partition, 1, mtx, mtx.shape[1] - ns)[:, ::-1][:, :ns[-1]]
    values = np.zeros((mtx.shape[0], len(ns)))
    acc = np.zeros((mtx.shape[0]))
    prev = 0
    for j, n in enumerate(ns):
        acc += partitioned[:, prev:n].sum(axis=1)
        values[:, j] = acc
        prev = n
    return values / sums[:, None]

@numba.jit(parallel=True)
def top_segment_proportions_sparse_csr(data, indptr, ns):
    ns = np.sort(ns)
    maxidx = ns[-1]
    sums = np.zeros((indptr.size - 1), dtype=data.dtype)
    values = np.zeros((indptr.size-1, len(ns)), dtype=np.float64)
    # Just to keep it simple, as a dense matrix
    partitioned = np.zeros((indptr.size-1, maxidx), dtype=data.dtype)
    for i in numba.prange(indptr.size - 1):
        start, end = indptr[i], indptr[i+1]
        sums[i] = np.sum(data[start:end])
        if end - start <= maxidx:
            partitioned[i, :end-start] = data[start:end]
        elif (end - start) > maxidx:
            partitioned[i, :] = - \
                (np.partition(-data[start:end], maxidx))[:maxidx]
    partitioned = np.apply_along_axis(
        np.partition, 1, partitioned, maxidx - ns)[:, ::-1][:, :ns[-1]]
    acc = np.zeros((indptr.size-1), dtype=data.dtype)
    prev = 0
    for j, n in enumerate(ns):
        acc += partitioned[:, prev:n].sum(axis=1)
        values[:, j] = acc
        prev = n
    return values / sums[:, None]
