"""
Template for Copy and Paste.
"""

import scanpy as sc

def myexample():
    """
    Template for preprocessing function. Use copy and paste.

    Returns
    -------
    adata : an AnnData object
    """

    # read a simple dictionary "ddata" that contains the data matrix, 
    # row_names and col_names from any path on your system
    path_to_data = 'data/myexample/'
    ddata = sc.read(path_to_data + 'myexample.csv')
    # other data reading examples
    #ddata = sc.read(path_to_data + 'myexample.txt')
    #ddata = sc.read(path_to_data + 'myexample.h5', sheet='mysheet')
    #ddata = sc.read(path_to_data + 'myexample.xlsx', sheet='mysheet')
    #ddata = sc.read(path_to_data + 'myexample.txt.gz')
    #ddata = sc.read(path_to_data + 'myexample.soft.gz')
    # if the first column does not store strings, rownames are not detected
    #  automatically, hence
    #ddata = sc.read(path_to_data + 'myexample.csv', first_column_names=True) 

    # from the simple dictionary, generate an AnnData object, which is similar
    # to R's ExpressionSet (Huber et al., Nat. Meth. 2015)
    # AnnData allows annotation of samples/cells and variables/genes via
    # the attributes "smp" and "var"
    adata = sc.AnnData(ddata)
    # adata = adata.transpose() # rows = samples/cells & columns = variables/genes
    
    # read some annotation from a file, now we want strings, and not a numerical
    # data matrix
    dgroups = sc.read(path_to_data + 'mygroups.csv', as_strings=True)
    adata.smp['groups'] = dgroups['X'][:, 0]

    # as with a dict, you can add arbitrary additional data to an data
    # for example, DPT needs a the expression vector of a root cell
    adata['xroot'] = adata.X[336]

    return adata


