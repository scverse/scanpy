
import numpy as np
from packaging import version
from scipy.sparse import csr_matrix
from sklearn import __version__ as sklearn_version
from sklearn.utils import check_random_state
from sklearn.utils.extmath import randomized_svd

# In older versions sklearn csr_mean_variance_axis0 always copies the whole sparse matrix
# if its dtype is float64
# https://github.com/scikit-learn/scikit-learn/blob/f0ab589f1541b1ca4570177d93fd7979613497e3/sklearn/utils/sparsefuncs_fast.pyx#L72

if version.parse(sklearn_version) < version.parse('0.20.0'):
    from sklearn.utils.sparsefuncs_fast import _csr_mean_variance_axis0
    mean_variance = lambda X: _csr_mean_variance_axis0(X.data, X.shape, X.indices)
else:
    from sklearn.utils.sparsefuncs_fast import csr_mean_variance_axis0 as mean_variance

#need to pass issparse check
class CentSparse(csr_matrix):
    def __init__(self, A, m=None, var=None):
        self.A = A
        self.m = m
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
    @property
    def dtype(self):
        return self.A.dtype
    @property
    def shape(self):
        return self.A.shape
    @property
    def T(self):
        return CentSparse(self.A.T, self.m.T)

class SparseDataPCA:
    def __init__(self, n_components=2, n_iter=5, random_state=None):
        self.n_components = n_components
        self.n_iter = n_iter
        self.random_state = random_state
    def fit_transform(self, X):
        k = self.n_components
        n_features = X.shape[1]
        r_s = check_random_state(self.random_state)

        m, var = mean_variance(X)
        m = np.asmatrix(m)
        var = var.sum()
        C = CentSparse(X, m)

        if k >= n_features:
            raise ValueError('n_components must be < n_features;'
                             ' got %d >= %d' % (k, n_features))
        U, Sigma, VT = randomized_svd(C, k, n_iter=self.n_iter, random_state=r_s)
        self.components_ = VT
        X_transformed = U * Sigma
        self.explained_variance_ = exp_var = np.var(X_transformed, axis=0)
        self.explained_variance_ratio_ = exp_var / var
        self.singular_values_ = Sigma
        return X_transformed
