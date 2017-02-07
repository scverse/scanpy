# Copyright 2016-2017 F. Alexander Wolf (http://falexwolf.de).
"""
Preprocess Data

Normalization and filtering functions that can be called directly or via
>>> ddata_or_X = preprocess(key, ddata_or_X)
or euqivalently
>>> ddata_or_X = pp(key, ddata_or_X)

Here, ddata_or_X is a data dictionary or data matrix and key is a string that
identifies the preprocessing function.
"""

def preprocess(key, ddata_or_X, *args, **kwargs):
    """
    Preprocess data with a set of available functions in this module.

    A function is selected based on the 'preprocess key' and
    is passed the arguments args and kwargs.

    Parameters
    ----------
    key : str
       String that identifies a preprocessing function.
    ddata_or_X : dict or np.ndarray
       Data dictionary containing at least a data matrix 'X'.

    Returns
    -------
    ddata_or_X : dict or np.ndarray
       Data dictionary that stores the preprocessed data matrix.
    """

    if isinstance(ddata_or_X, dict):
        isddata = True
    else:
        isddata = False

    try:
        if isddata:
            from . import advanced
            ppfunc = getattr(advanced, key)
        else:
            from . import simple
            ppfunc = getattr(simple, key)
    except AttributeError:
            msg = ('Do not know how to run preprocessing key "' + key
                   + '" for argument "'
                   + ('ddata' if isddata else 'X')
                   + '".\nEither define a function ' + key + '() '
                   'in scanpy/preprocess/' 
                   + ('advanced' if isddata else 'simple') 
                   + '.\nOr use one of the available functions\n'
                   + ppkeys_str())
            from sys import exit
            exit(msg)
    return ppfunc(ddata_or_X, *args, **kwargs)

pp = preprocess
"""
Same as preprocess().
"""

def ppkeys_str():
    from . import advanced, simple
    str = '... that take ddata as argument (advanced)'
    for k in sorted(dir(advanced)):
        if not k.startswith('_'):
            str += '\n    ' + k
    str += '\n... that take X as argument (simple)'
    for k in sorted(dir(simple)):
        if not k.startswith('_'):
            str += '\n    ' + k
    return str
