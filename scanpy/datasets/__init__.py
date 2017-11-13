# Author: F. Alex Wolf (http://falexwolf.de)
"""Init runs, manage examples.
"""


from .builtin import *
from ..utils import sanitize_anndata


def init_run(run_name, suffix='', recompute=True, reread=False,
             return_module=False):
    """Read and preprocess data based on a "run file".

    Filenames of the form "runs_whatevername.py", "scanpy_whatevername.py" and
    "preprocessing_whatevername.py" in the current working directory are
    automatically considered as run files.

    In addition, there are builtin examples defined in the run file
    https://github.com/theislab/scanpy/tree/master/scanpy/examples/builtin.py

    Parameters
    ----------
    run_name : str
        Key for looking up an example-preprocessing function.
    suffix : str, optional (default: '')
        Set suffix to be appended to `run_name` in naming output files.
    recompute : bool, optional (default: True)
        Recompute preprocessing.
    reread : bool, optional (default: False)
        Reread the original data file (often a text file, much slower) instead
        of the hdf5 file.
    return_module : bool, optional (default: False)
        Return example module.

    Returns
    -------
    adata : AnnData
        Annotated data matrix, optionally with metadata such as
        adata.uns['xroot'] : np.ndarray or int
            Expression vector or index of root cell for DPT analysis.
    Additionally, if return_module == True:
    exmodule : dict, optional
        Example module.
    """
    import os, sys
    from .. import readwrite
    from .. import settings as sett
    from .. import logging as logg
    sett._run_basename = run_name
    sett._run_suffix = suffix
    sett.run_name = sett._run_basename + sett._run_suffix
    adata_file = readwrite.get_filename_from_key(sett.run_name)
    adata_file_exists = os.path.exists(adata_file)
    # find the runfile with preprocessing functions etc.
    loop_over_filenames = [filename for filename in os.listdir('.')
                           if (filename.startswith('runs')
                               or filename.startswith('preprocessing')
                               or filename.startswith('scanpy'))
                           and filename.endswith('.py')]
    if len(loop_over_filenames) == 0:
        logg.m('did not find user examples, to provide some,\n'
               '    generate a file "runs_whatevername.py" in your working directory,\n'
               '    analogous to https://github.com/theislab/scanpy/blob/master/scanpy/examples/builtin.py',
               v='hint')
    not_found = True
    sys.path.insert(0, '.')
    for filename in loop_over_filenames:
        exmodule = __import__(filename.replace('.py', ''))
        try:
            exfunc = getattr(exmodule, run_name)
            not_found = False
        except AttributeError:
            pass
    if not_found:
        try:
            exfunc = getattr(builtin, run_name)
            exmodule = builtin
        except AttributeError:
            import types
            sys.exit('Do not know how to run example "{}".\nEither define a function {}() '
                     'that returns an AnnData object in "./runfile_whatevername.py".\n'
                     'Or, use one of the builtin examples: {}'
                     .format(run_name, run_name,
                             [a for a in dir(builtin)
                              if isinstance(getattr(builtin, a), types.FunctionType)]))
    if not adata_file_exists or recompute or reread:
        logg.m('reading and preprocessing data')
        # run the function
        adata = exfunc()
        # add run_name to adata
        logg.m('... X has shape n_samples × n_variables = {} × {}'
               .format(adata.X.shape[0], adata.X.shape[1]))
        # do sanity checks on data dictionary
        adata = sanitize_anndata(adata, verbosity=1)
        # write the prepocessed data
        readwrite.write(sett.run_name, adata)
    else:
        adata = readwrite.read(sett.run_name)

    if return_module:
        return adata, exmodule
    else:
        return adata


# -------------------------------------------------------------------------------
# Reading and writing with sett.run_name
# -------------------------------------------------------------------------------


def read_run(run_name=None, suffix=''):
    """Read run and init sett.run_name if provided.
    """
    from .. import settings as sett
    if run_name is None: run_name = sett.run_name
    if suffix == '': suffix = sett._run_suffix
    sett._run_basename = run_name
    sett._run_suffix = suffix
    sett.run_name = sett._run_basename + sett._run_suffix
    return init_run(run_name, suffix=suffix, recompute=False)


def write_run(data, ext=None):
    """Write run.

    ext : str or None (default: None)
        File extension from wich to infer file format.
    """
    from .. import readwrite
    from .. import settings as sett
    readwrite.write(sett.run_name, data, ext=ext)
