from __future__ import annotations

from importlib import import_module
from types import ModuleType
from typing import Literal, Union
from collections.abc import Generator

from sklearn.base import TransformerMixin

from ._common import mappings, select_backend


_Eponymous = Literal['annoy', 'faiss', 'nmslib', 'pynndescent']
_Backend = Union[_Eponymous, Literal['sklearn']]
_Algorithm = Union[_Eponymous, Literal['ball_tree', 'kd_tree', 'brute']]

_backends: dict[_Backend, set[_Algorithm]] = {
    'sklearn': {'ball_tree', 'kd_tree', 'brute'},
    'annoy': {'annoy'},
    'faiss': {'faiss'},
    'nmslib': {'nmslib'},
    'pynndescent': {'pynndescent'},
}

BACKENDS, ALGORITHMS = mappings(_backends)


def import_backend(backend: _Backend) -> ModuleType:
    return import_module(f'sklearn_ann.kneighbors.{backend}')


def get_transformer(
    algorithm: _Algorithm, backend: _Backend | None = None
) -> type[TransformerMixin]:
    actual_backend: _Backend = select_backend(ALGORITHMS, algorithm, backend)
    mod = import_backend(actual_backend)
    transformer_cls: type[TransformerMixin] = next(
        cls
        for name, cls in vars(mod).items()
        if name.lower().startswith(algorithm.replace('_', ''))
    )
    return transformer_cls


def available_backends() -> Generator[_Backend, None, None]:
    for backend in BACKENDS:
        try:
            import_backend(backend)
        except ImportError:
            pass
        else:
            yield backend
