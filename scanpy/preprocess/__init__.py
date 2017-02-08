# Copyright 2016-2017 F. Alexander Wolf (http://falexwolf.de).
"""
Preprocess Data
"""

from .simple import *
from .advanced import *

__all__ = [
    # from simple
    'filter_cells',
    'filter_genes_cv',
    'filter_genes_fano',
    'log',
    'pca',
    'row_norm',
    'zscore',
    # from advanced
    'subsample',
    'weinreb16'
]

def show():
    from . import advanced, simple
    str = 'take ddata as argument (advanced)'
    for k in sorted(dir(advanced)):
        if not k.startswith('_'):
            str += '\n    ' + k
    str += '\ntake X as argument (simple)'
    for k in sorted(dir(simple)):
        if not k.startswith('_'):
            str += '\n    ' + k
    return str
