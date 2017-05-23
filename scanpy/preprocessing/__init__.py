# Author: F. Alex Wolf (http://falexwolf.de)
"""Preprocessing Functions

Simple functions and whole recipes.
"""

from .simple import *
from .recipes import *


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
