import numpy as np
from anndata import AnnData
from scanpy.api import Neighbors

def test_compute_similarities():
    # the result
    similarities = [
        [0.0, 0.5146393179893494, 0.0, 0.36445462703704834],
        [0.5146393179893494, 0.0, 0.3581143319606781, 0.2239987552165985],
        [0.0, 0.3581143319606781, 0.0, 0.5245543718338013],
        [0.36445462703704834, 0.2239987552165985, 0.5245543718338013, 0.0]]
    # the input data
    X = np.array([[1, 0], [3, 0], [5, 6], [0, 4]])
    adata = AnnData(X)
    neighbors = Neighbors(adata)
    neighbors.compute_distances()
    neighbors.compute_similarities()
    np.allclose(
        neighbors.similarities.toarray(), similarities, np.finfo(np.float32).eps)
