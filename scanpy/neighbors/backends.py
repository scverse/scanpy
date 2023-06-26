from __future__ import annotations

from typing import Literal
from collections.abc import Mapping
from sklearn.neighbors import KNeighborsTransformer


_Backend = Literal['annoy', 'faiss', 'nmslib', 'pynndescent', 'sklearn']

ALGORITHMS: Mapping[_Backend, frozenset[str]] = {
    'sklearn': frozenset({'ball_tree', 'kd_tree', 'brute'}),
    'annoy': frozenset({'annoy'}),
    'faiss': frozenset({'faiss'}),
    'nmslib': frozenset({'nmslib'}),
    'pynndescent': frozenset({'pynndescent'}),
}


def get_transformer(backend: _Backend) -> KNeighborsTransformer:
    pass  # from sklearn_ann.kneighbors.sklearn import ...


def available_backends():
    pass
