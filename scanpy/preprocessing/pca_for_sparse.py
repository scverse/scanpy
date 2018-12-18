from scipy.sparse import csr_matrix
from sklearn.utils.sparsefuncs_fast import csr_mean_variance_axis0
import numpy as np

#need to pass issparse check
class CentSparse(csr_matrix):
    def __init__(self, A, m=None, var=None):
        self.A = A
        if m is not None:
            self.m = m
            self.var = var
        else:
            self.m, self.var = csr_mean_variance_axis0(A) #needs sklearn version 0.20.1 to be efficient
        self.c_r = self.m.shape[0] == 1
    def __mul__(self, B):
        AB = self.A * B
        mB = np.dot(self.m, B) if self.c_r else self.m * B.sum(0)[None, :]
        return AB - mB
    def __rmul__(self, B):
        BA = B * self.A
        Bm = B.sum(1)[:, None] * self.m if self.c_r else np.dot(B, self.m)
        return BA - Bm
    def toarray(self):
        return self.A - self.m
    def astype(self, dtype):
        if dtype == self.dtype:
            return self
        return CentSparse(self.A.astype(dtype), self.m.astype(dtype))
    @property
    def dtype(self):
        return self.A.dtype
    @property
    def shape(self):
        return self.A.shape
    @property
    def T(self):
        return CentSparse(self.A.T, self.m.T)
    @property
    def data(self):
        return self.A.data
    @property
    def indices(self):
        return self.A.indices
