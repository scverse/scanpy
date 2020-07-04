"""This file contains some common fixtures for use in tests.

This is kept seperate from the helpers file because it relies on pytest.
"""
import pytest
import numpy as np
from scipy import sparse

from anndata.tests.helpers import asarray


@pytest.fixture(
    params=[sparse.csr_matrix, sparse.csc_matrix, asarray],
    ids=["scipy-csr", "scipy-csc", "np-ndarray"],
)
def array_type(request):
    """Function which converts passed array to one of the common array types."""
    return request.param


@pytest.fixture(params=[np.float64, np.float32])
def float_dtype(request):
    return request.param
