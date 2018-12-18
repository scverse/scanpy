from scipy.sparse import csr_matrix
import numpy as np

#need to pass issparse check
class CentSparse(csr_matrix):
    def __init__(self, A, m=None):
        self.A = A
        self.m = A.mean(0) if m is None else m
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
