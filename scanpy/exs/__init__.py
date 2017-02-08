# Copyright 2016-2017 F. Alexander Wolf (http://falexwolf.de).
"""
Example Data and Example Use Cases
"""

from . import builtin, user
from .. import utils
from .. import readwrite
from .. import settings as sett

def exdata(format='plain'):
    """Show available example data.
    """
    if format == 'plain':
        s = utils.pretty_dict_string(dexdata())
    elif format == 'markdown':
        s = utils.markdown_dict_string(builtin.dexdata)
    print(s)

def examples():
    """Show available example use cases.
    """
    s = utils.pretty_dict_string(dexamples())
    print(s)

def dexdata(): 
    """Example data.
    """
    all_dex = utils.merge_dicts(builtin.dexdata, user.dexdata)
    try:
        # additional possibility to add example module
        from . import user_private
        all_dex = utils.merge_dicts(all_dex, user_private.dexdata) 
    except ImportError:
        pass
    return all_dex

def dexamples():
    """Example use cases.
    """
    builtin_dex = utils.fill_in_datakeys(builtin.dexamples, builtin.dexdata)
    user_dex = utils.fill_in_datakeys(user.dexamples, user.dexdata)
    all_dex = utils.merge_dicts(builtin_dex, user_dex) 
    try:
        # additional possibility to add example module
        from . import user_private
        user_private_dex = utils.fill_in_datakeys(user_private.dexamples, 
                                                  user_private.dexdata)
        all_dex = utils.merge_dicts(all_dex, user_private_dex) 
    except ImportError:
        pass
    return all_dex

def example(exkey, return_module=False):
    """
    Read and preprocess data for predefined example.

    Parameters
    ----------
    exkey : str
        Key for the example dictionary in _examples.
    return_module : bool, optional
        Return example module.

    Returns
    -------
    adata : dict containing
        X : np.ndarray
            Data array for further processing, columns correspond to genes,
            rows correspond to samples.
        row_names : np.ndarray
            Array storing the experimental labels of samples.
        col_names : np.ndarray
            Array storing the names of genes.

        There might be further entries such as
        groupnames_n : np.ndarray
            Array of shape (number of samples) that indicates groups.
        xroot : np.ndarray or int
            Expression vector or index of root cell for DPT analysis.

    Returns additionally, if return_module == True:
    exmodule : dict, optional
        Example module.
    """
    try:
        exfunc = getattr(user, exkey)
        exmodule = user
    except AttributeError:
        try:
            # additional possibility to add example module
            from . import user_private
            exfunc = getattr(user_private, exkey)
            exmodule = user_private
        except (ImportError, AttributeError):
            try:
                exfunc = getattr(builtin, exkey)
                exmodule = builtin
            except AttributeError:
                msg = ('Do not know how to run example "' + exkey +
                       '".\nEither define a function ' + exkey + '() '
                       'in scanpy/exs/user.py that returns a data dictionary.\n'
                       'Or use one of the available examples:'
                       + exkeys_str())
                from sys import exit
                exit(msg)

    from os.path import exists
    exfile = readwrite.get_filename_from_key(sett.basekey)
    if (not exists(exfile)
        or sett.recompute in ['read', 'pp']):
        # run the function
        adata = exfunc()
        # add exkey to adata
        adata['exkey'] = exkey
        sett.m(0, 'X has shape nr_samples x nr_variables =', 
               adata.X.shape[0], 'x', adata.X.shape[1])
        # do sanity checks on data dictionary
        adata = check_adata(adata)
        readwrite.write(sett.basekey, adata)
        sett.m(0, 'wrote preprocessed data to', exfile)
    else:
        adata = readwrite.read(sett.basekey)

    if return_module:
        return adata, exmodule
    else:
        return adata

#-------------------------------------------------------------------------------
# Checking of data dictionary
# - Might be replaced with a data class.
#-------------------------------------------------------------------------------


ignore_groups = ['N/A', 'dontknow', 'no_gate']
# howtos
howto_specify_subgroups = '''sample annotation in adata only consists of sample names
--> you can provide additional annotation by setting, for example,
    adata.smp['groups'] = ['A', 'B', 'A', ... ]
    adata.smp['time'] = [0.1, 0.2, 0.7, ... ]'''

def check_adata(adata):
    """
    Do sanity checks on adata object.

    Checks whether adata contains annotation.
    """
    import numpy as np
    import sys
    if len(adata.smp_keys()) == 0:
        sett.m(0, howto_specify_subgroups)
    else:
        for k in adata.smp_keys():
            # ordered unique categories
            if not k + '_names' in adata:
                adata[k + '_names'] = np.unique(adata.smp[k])
            adata[k + '_names'] = np.setdiff1d(adata[k + '_names'],
                                               np.array(ignore_groups))
            # output 
            sett.m(0,'sample annotation', k, 'with', adata[k + '_names'])
            # indices for each category
            if not k + '_ids' in adata:
                adata[k + '_ids'] = np.arange(len(adata[k + '_names']), dtype=int)
            # masks for each category
            if not k + '_masks' in adata:
                masks = []
                for name in adata[k + '_names']:
                    masks.append(name == adata.smp[k])
                adata[k + '_masks'] = np.array(masks)
    return adata

howto_specify_subgroups_ddata = '''no key "row" in adata dictionary found
--> you might provide a dict of lists of n subgroup names (strings or ints) with
    number of samples as follows
    {'group1': ['A', 'B', 'A', ... ], 'group2': ['c', 'a', ...]}
'''

def check_ddata(ddata):
    """
    Do sanity checks on ddata dictionary.

    Checks whether ddata conains categorical row metadata 'row'. 

    If yes, for each class of categories in 'row' associate an 'order',
    indices 'ids', colors and masks.
    """
    import numpy as np
    import sys
    if not 'row' in ddata:
        sett.m(0, howto_specify_subgroups_ddata)
    else:
        if not isinstance(ddata['row'], dict):
            msg = 'row must be a dictionary! {\'cat1\': [...], }'
            sys.exit(msg)
        for k in ddata['row']:
            # transform to np.ndarray
            try:
                ddata['row'][k] = np.array(ddata['row'][k], dtype=int)
            except:
                ddata['row'][k] = np.array(ddata['row'][k], dtype=str)
            # ordered unique categories
            if not k + '_names' in ddata:
                ddata[k + '_names'] = np.unique(ddata['row'][k])
            # output 
            sett.m(0,'read sample annotation', k, 'with', ddata[k + '_names'])
            # indices for each category
            if not k + '_ids' in ddata:
                ddata[k + '_ids'] = np.arange(len(ddata[k + '_names']), dtype=int)
            # masks for each category
            if not k + '_masks' in ddata:
                masks = []
                for name in ddata[k + '_names']:
                    masks.append(name == ddata['row'][k])
                ddata[k + '_masks'] = np.array(masks)
    return ddata

def exkeys_str():
    str = '\n'
    for k in sorted(dexamples().keys()):
        str += '    ' + k + '\n'
    return str
