import cython
import numpy as np
cimport numpy as np
from libc.math cimport sqrt

# The following did not lead to a speed-up.
# @cython.boundscheck(False)
# def get_Ddiff_row(int i, np.ndarray[np.float32_t, ndim=1] evals,
#               np.ndarray[np.float32_t, ndim=2] rbasis,
#               np.ndarray[np.float32_t, ndim=2] lbasis):
#     cdef np.ndarray[np.float32_t, ndim=1] m_i = get_M_row(i, evals, rbasis, lbasis)
#     cdef np.ndarray[np.float32_t, ndim=1] d_i = np.zeros([rbasis.shape[0]], dtype=np.float32)
#     cdef int j
#     for j in range(rbasis.shape[0]):
#         m_j = get_M_row(j, evals, rbasis, lbasis)
#         d_i[j] = c_dist(m_i, m_j)
#     return d_i

@cython.boundscheck(False)
def get_M_row(int i, np.ndarray[np.float32_t, ndim=1] evals,
              np.ndarray[np.float32_t, ndim=2] rbasis,
              np.ndarray[np.float32_t, ndim=2] lbasis):
    cdef np.ndarray[np.float32_t, ndim=1] m_i = np.zeros([rbasis.shape[0]], dtype=np.float32)
    cdef np.float32_t factor
    cdef int l
    cdef int j
    for l in range(1, evals.shape[0]):
        factor = evals[l]/(1-evals[l])
        for j in range(lbasis.shape[0]):
            m_i[j] += factor * rbasis[i, l] * lbasis[j, l]
    m_i += rbasis[i, 0] * lbasis[:, 0]
    return m_i

@cython.boundscheck(False)
def c_dist(np.ndarray[np.float32_t, ndim=1] vec1,
           np.ndarray[np.float32_t, ndim=1] vec2):
    cdef np.float32_t dist = 0
    cdef np.float32_t tmp
    cdef int i
    for i in range(vec1.size):
        tmp = vec1[i] - vec2[i]
        dist +=  tmp * tmp
    return sqrt(dist)
