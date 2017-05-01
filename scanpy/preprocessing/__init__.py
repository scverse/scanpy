# Copyright 2016-2017 F. Alexander Wolf (http://falexwolf.de).
"""
Preprocess Data
"""

from .simple import *
from .recipes import *

__all__ = [
    # from simple
    'filter_cells',
    'gene_cv_filter',
    'filter_genes_fano',
    'log',
    'pca',
    'row_norm',
    'zscore',
    # from recipes
    'subsample',
    'weinreb16'
]

def overview():
    from . import recipes, simple
    str = 'take adata as argument (recipes)'
    for k in sorted(dir(recipes)):
        if not k.startswith('_'):
            str += '\n    ' + k
    str += '\ntake X as argument (simple)'
    for k in sorted(dir(simple)):
        if not k.startswith('_'):
            str += '\n    ' + k
    return str
