"""
Template for Copy and Paste.
"""

import numpy as np
import scanpy as sc

def myexample():
    """
    Template for preprocessing function. Use copy and paste.

    Returns
    -------
    adata : AnnData
        Stores data matrix and sample and variable annotations as well
        as an arbitrary amount of unstructured annotation. For the latter
        it behaves like a Python dictionary.
    """
    # Generate an AnnData object, which is similar
    # to R's ExpressionSet (Huber et al., Nat. Meth. 2015)
    # AnnData allows annotation of samples/cells and variables/genes via
    # the attributes "smp" and "var"
    path_to_data = 'data/myexample/'
    adata = sc.read(path_to_data + 'myexample.csv')
    # other data reading examples
    #adata = sc.read(path_to_data + 'myexample.txt')
    #adata = sc.read(path_to_data + 'myexample.h5', sheet='mysheet')
    #adata = sc.read(path_to_data + 'myexample.xlsx', sheet='mysheet')
    #adata = sc.read(path_to_data + 'myexample.txt.gz')
    #adata = sc.read(path_to_data + 'myexample.soft.gz')
    # if the first column does not store strings, rownames are not detected
    # automatically, hence
    #adata = sc.read(path_to_data + 'myexample.csv', first_column_names=True)

    # transpose if needed to match the convention that rows store samples/cells
    # and columns variables/genes
    # adata = adata.transpose() # rows = samples/cells & columns = variables/genes

    # read some annotation from a file, now we want strings, and not a numerical
    # data matrix, the following reads from the first column of the file
    groups = np.genfromtext(path_to_data + 'mygroups.csv', dtype=str)
    adata.smp['groups'] = groups[:, 0]
    # or alternatively, when you want to be smart about row and column annotaton
    # dgroups = sc.read(path_to_data + 'mygroups.csv', as_strings=True, return_dict=True)
    # adata.smp['groups'] = dgroups['X'][:, 0]

    # as with a dict, you can add arbitrary additional data to an data
    # for example, DPT needs a the expression vector of a root cell
    adata['xroot'] = adata.X[336]
    return adata


