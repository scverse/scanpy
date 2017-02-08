"""
Example Data and Use Cases - User Examples

Add example data and use cases while being assured that no conflicts with
scanpy/exs/builtin.py in the master branch on GitHub arise. The latter
might be updated by other contributors.

If you want to update the builtin Scanpy examples on GitHub, copy and paste the
example from here to there and make a pull request.

Attributes
----------
Functions read data and do preprocessing, furthermore

dexdata : dict
    Stores information about example data.
dexamples : dict 
    Stores information about example use cases. The keys in this dictionary also
    exist as a function attribute of this module that performs data reading and
    preprocessing. 
"""

# this is necessary to import scanpy from within package
from __future__ import absolute_import
# scientific modules
import numpy as np
# scanpy
import scanpy as sc
from .. import utils
from .. import settings as sett

#--------------------------------------------------------------------------------
# The 'dexdata dictionary' stores information about example data.
# - is optional, can stay empty
#--------------------------------------------------------------------------------

dexdata = {
}

#--------------------------------------------------------------------------------
# The 'example dictionary' provides information about tool parameters 
# that deviate from default parameters.
# - is optional, can stay empty
#--------------------------------------------------------------------------------

dexamples = {
}

#--------------------------------------------------------------------------------
# One function per example that reads, annotates and preprocesses data
# - one function 'exkey()' per 'exkey'
#--------------------------------------------------------------------------------

def myexample():

    # get help
    # help(sc.read)

    # read data from any path on your system
    path_to_data = 'data/myexample/'
    ddata = sc.read(path_to_data + 'myexample.csv')
    adata = AnnData(ddata)

    # other data reading examples
    # adata = sc.read(path_to_data + 'myexample.csv', first_column_names=True)
    # adata = sc.read(path_to_data + 'myexample.h5', sheet='countmatrix')
    # adata = sc.read(path_to_data + 'myexample.xlsx', sheet='countmatrix')
    # adata = sc.read(path_to_data + 'myexample.txt', sheet='countmatrix')
    # adata = sc.read(path_to_data + 'myexample.txt.gz', sheet='countmatrix')
    # adata = sc.read(path_to_data + 'myexample.soft.gz', sheet='countmatrix')

    # in adata.X, rows should correspond to samples, columns to genes
    # to match this convention, transpose your data if necessary
    # adata = utils.transpose_adata(adata)

    # get groupnames (as strings)
    dgroups = sc.read(path_to_data + 'mygroups.csv', as_strings=True)
    adata['groupnames_n'] = dgroups['X'][:, 0]

    # specify root cell
    adata['xroot'] = adata.X[336]

    return adata

#--------------------------------------------------------------------------------
# Optional functions for Raw Data, Annotation, Postprocessing, respectively
#--------------------------------------------------------------------------------

