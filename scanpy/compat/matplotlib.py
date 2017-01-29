from __future__ import absolute_import  # else this module tries to import itself

import warnings

__all__ = ['pyplot']

with warnings.catch_warnings():
    # import of pyplot raises warning only in Python 2
    warnings.simplefilter('ignore')
    from matplotlib import pyplot
