# Copyright 2016-2017 F. Alexander Wolf (http://falexwolf.de).
"""
Example Data and Example Use Cases
"""

from . import builtin
from .. import utils
from .. import readwrite
from .. import settings as sett

def example(exkey, subsample=1, return_module=False):
    """
    Read and preprocess data for predefined example.

    Parameters
    ----------
    exkey : str
        Key for the example dictionary in _examples.
    subsample : int
        Subsample to a fraction of 1/subsample of the data.
    return_module : bool, optional
        Return example module.

    Returns
    -------
    adata : AnnData
        Annotated data matrix, optionally with metadata such as
        adata['xroot'] : np.ndarray or int
            Expression vector or index of root cell for DPT analysis.

    Returns additionally, if return_module == True:
    exmodule : dict, optional
        Example module.
    """
    try:
        try:
            from sys import path
            path.insert(0, '.')
            import user_exs
        except ImportError:
            sett.m(0, '--> did not find user examples, to provide some,\n'
                   '    generate file user_exs.py in your working directory')
        exfunc = getattr(user_exs, exkey)
        exmodule = user_exs
    except (UnboundLocalError, AttributeError):
        try:
            # additional possibility to add example module
            from . import builtin_private
            exfunc = getattr(builtin_private, exkey)
            exmodule = builtin_private
        except (ImportError, AttributeError):
            try:
                exfunc = getattr(builtin, exkey)
                exmodule = builtin
            except AttributeError:
                msg = ('Do not know how to run example "' + exkey +
                       '".\nEither define a function ' + exkey + '() '
                       'in ./user_exs.py that returns an AnnData object.\n'
                       'Or, use one of the builtin examples:'
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
        # adata['exkey'] = exkey
        sett.m(0, 'X has shape n_samples x n_variables =', 
               adata.X.shape[0], 'x', adata.X.shape[1])
        # do sanity checks on data dictionary
        adata = check_adata(adata)

        # subsampling
        if subsample != 1:
            from ..preprocess import subsample as subsample_function
            adata = subsample_function(adata, subsample)

        readwrite.write(sett.basekey, adata)
        sett.m(0, 'wrote preprocessed data to', exfile)
    else:
        adata = readwrite.read(sett.basekey)

    if return_module:
        return adata, exmodule
    else:
        return adata

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
    all_dex = utils.merge_dicts(builtin.dexdata, {})
    try:
        # additional possibility to add example module
        from . import builtin_private
        all_dex = utils.merge_dicts(all_dex, builtin_private.dexdata) 
    except ImportError:
        pass
    return all_dex

def dexamples():
    """Example use cases.
    """
    builtin_dex = utils.fill_in_datakeys(builtin.dexamples, builtin.dexdata)
    all_dex = utils.merge_dicts(builtin_dex, {}) 
    try:
        # additional possibility to add example module
        from . import builtin_private
        builtin_private_dex = utils.fill_in_datakeys(builtin_private.dexamples, 
                                                  builtin_private.dexdata)
        all_dex = utils.merge_dicts(all_dex, builtin_private_dex) 
    except ImportError:
        pass
    return all_dex

#-------------------------------------------------------------------------------
# Checks of AnnData object
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
    if 'tools' not in adata:
        adata['tools'] = np.array([], dtype=str)
    if len(adata.smp_keys()) == 0:
        sett.m(0, howto_specify_subgroups)
    else:
        for smp in adata.smp_keys():
            # ordered unique categories
            if not smp + '_names' in adata:
                adata[smp + '_names'] = np.unique(adata.smp[smp])
                try:
                    from natsort import natsorted
                    adata[smp + '_names'] = np.array(natsorted(adata[smp + '_names']))
                except:
                    pass
                if adata[smp + '_names'].dtype.char == 'U':
                    adata[smp + '_names'] = np.setdiff1d(adata[smp + '_names'],
                                                       np.array(ignore_groups))
            # output
            ann_info = adata[smp + '_names']
            if len(adata[smp + '_names']) > 20:
                ann_info = (str(adata[smp + '_names'][0:3]).replace(']','') 
                            + ' ...' + str(adata[smp + '_names'][-2:]).replace('[',''))
            sett.m(0,'sample annotation','"' + smp + '"', 'with',  ann_info)
    return adata

def exkeys_str():
    str = ''
    for k in sorted(dexamples().keys()):
        str += '\n    ' + k
    return str
