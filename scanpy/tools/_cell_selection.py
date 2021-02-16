"""
Returns a selection of cells simultaneously expressing given genes. 
"""
from typing import Optional, Sequence, Union

import pandas as pd
from anndata import AnnData

from .. import logging as logg
from .._compat import Literal

def cell_selection(
    adata: AnnData,
    var_names: Sequence[str] = None,
    denoise: bool = True,
    t: Union[Literal['auto'], int] = 'auto',
    knn: int = 10,
    key_added: Optional[str] = None,
    inplace: bool = True,
    **kwargs,
) -> Optional[Anndata]:
    """\
    Returns a selection of cells simultaneously expressing given genes. 

    The selected cells must have a nonnegative expression of all the supplied genes.

    By default the data is denoised by expression imputation with MAGIC. 

    Parameters
    ----------
    adata
        Annotated data matrix
    var_names
        List of var_names to use for cell selection.
    denoise
        Whether to use imputation method to denoise the data (recommended).
        Requires the MAGIC package.
    t
        Parameter for MAGIC denoising. 
        Power to which the diffusion operator is powered. 
        See :func:`scanpy.external.pp.magic` for more information.
    knn
        Parameter for MAGIC denoising. 
        Number of nearest neighbors on which to build kernel.
        See :func:`scanpy.external.pp.magic` for more information.
    key_added
        By default, the selection information is added to
        `.obs[f'select_{"_".join(var_names)}']`.
        Notice that the `var_names` labels are added to the cell selection name.
    inplace
        If `True`, adds selection information to `adata.obs[key_added]`,
        else this function returns the information.
    kwargs
        Additional arguments to `magic.MAGIC`.

    Returns
    -------
    A pandas dataframe with the selection of cells expressing given genes if `inplace=False`.
    For `inplace=True` `adata.obs` is updated with an additional field
    specified by the `key_added` parameter (default = 'select_{"_".join(var_names)').


    Examples
    --------
    >>> import scanpy as sc
    >>> adata = sc.datasets.pbmc68k_reduced()
    >>> gene_list = ['CD14', 'LYZ']
    >>> sc.tl.cell_selection(adata, var_names=gene_list)
    >>> sc.pl.umap(adata, color='select_CD14_LYZ', groups=[True])
    """

    for var in var_names:
            if var not in adata.var_names:
                raise ValueError(
                    f'Given value: {var} is not a valid name from {adata.var_names}'
                )


    if denoise:  
        try:
            from .. external.pp import magic as magic
        except ImportError:
            raise ImportError(
            'Please install magic package via `pip install --user '
            'git+git://github.com/KrishnaswamyLab/MAGIC.git#subdirectory=python`'
        )     
        adata_magic = magic(adata, name_list=var_names, t=t, knn=knn)

        val = adata_magic[:,'{}'.format(var_names[0])].X > 0
        for i in range(len(var_names[1:])):
            val = val & (adata_magic[:,'{}'.format(var_names[i+1])].X > 0)

        val_list = [item for sublist in val.tolist() for item in sublist]
        dat = pd.Categorical(val_list, categories=[True, False])
    else:
        val = adata.raw[:,'{}'.format(var_names[0]))].X.todense() > 0
        for i in range(len(var_names[1:])):
            val = val & (adata.raw[:,'{}'.format(var_names[i+1]))].X.todense() > 0)

        val_list = [item for sublist in val.tolist() for item in sublist]
        dat = pd.Categorical(val_list, categories=[True, False])

    if inplace:
        if key_added is None:
            key_added = f'select_{"_".join(var_names)}'
        logg.info(f'Storing cell selection info using `.obs[{key_added!r}]`')
        adata.obs[key_added] = dat
    else:
        return dat
