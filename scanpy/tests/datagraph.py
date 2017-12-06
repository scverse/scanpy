import numpy as np
import pandas as pd
from anndata import AnnData

from scipy.sparse import csr_matrix
from scanpy.api import DataGraph

def test_compute_transition_matrix():
    # the result
    Ktilde_result = [
        [0.0, 0.5146393179893494, 0.0, 0.36445462703704834],
        [0.5146393179893494, 0.0, 0.3581143319606781, 0.2239987552165985],
        [0.0, 0.3581143319606781, 0.0, 0.5245543718338013],
        [0.36445462703704834, 0.2239987552165985, 0.5245543718338013, 0.0]]
    # the input data
    X = np.array([[1, 0], [3, 0], [5, 6], [0, 4]])
    adata = AnnData(X)
    g = DataGraph(adata, n_jobs=1)
    g.compute_transition_matrix()
    np.allclose(g.Ktilde.toarray(), Ktilde_result, np.finfo(np.float32).eps)
