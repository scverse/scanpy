from scipy.sparse import csr_matrix

class CentSparse(csr_matrix):
    def __init__(self, A):
        self.A = A
        self.m = A.mean(0).A1
    def __mul__(self, B):
        return self.A * B - m * B
    def __rmul__(self, B):
        return B * self.A - B * m
