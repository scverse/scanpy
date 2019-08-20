"""\
Annotates gene expression (cell data) with cell types.
"""
import warnings

from anndata import AnnData
import pandas as pd


def annotator(
    adata: AnnData,
    markers: pd.DataFrame,
    num_genes: int = None,
    return_nonzero_annotations: bool = True,
    p_threshold: float = 0.05,
    p_value_fun: str = "binom",
    z_threshold: float = 1.0,
    scoring: str = "exp_ratio",
    normalize: bool = False,
):
    """\
    Annotator marks the data with cell type annotations based on marker genes.

    Over-expressed genes are selected with the  Mann-Whitney U tests and cell
    types are assigned with the hypergeometric test. This function first selects
    genes from gene expression data with the Mann Whitney U test, then annotate
    them with the hypergeometric test, and finally filter out cell types that
    have zero scores for all cells. The results are scores that tell how
    probable is each cell type for each cell.

    Parameters
    ----------
    adata : AnnData
        Tabular data with gene expressions.
    markers : pd.DataFrame
        The data-frame with marker genes and cell types. Data-frame has two
        columns **Gene** and **Cell Type** first holds gene names or ID and
        second cell type for this gene. Gene names must be written in the same
        format than genes in `adata`.
    num_genes : int, optional (default=None)
        The number of genes that the organism has.
    return_nonzero_annotations : bool, optional (default=True)
        If true return scores only for cell types that have no zero scores.
    p_threshold : float
        A threshold for accepting the annotations. Annotations that have FDR
        value bellow this threshold are used.
    p_value_fun : str, optional (defaults: binom)
        A function that calculates a p-value. It can be either
        `binom` that uses binom.sf or
        `hypergeom` that uses hypergeom.sf.
    z_threshold: float
        The threshold for selecting the gene from gene expression data.
        For each cell the attributes with z-value above this value are selected.
    scoring : str, optional (default = "exp_ratio")
        Scoring method for cell type scores. Available scores are:
        - exp_ratio - the proportion of genes typical for a cell type expressed
        in the cell
        - sum_of_expressed_markers - the sum of expressions of genes typical for
        a cell type
        - log_fdr - negative of the logarithm of an false discovery rate (FDR)
        value
        - log_p_value - negative of the logarithm of an p-value
    normalize : bool, optional (default = False)
        If this parameter is True data will be normalized during the
        a process with a log CPM normalization.
       That method works correctly data needs to be normalized.
        Set this `normalize` on True if your data are not normalized already.

    Returns
    -------
    pd.DataFrame
        Cell type for each cell for each cell. The result is a sore matrix that
        tells how probable is each cell type for each cell. Columns are cell
        types and rows are cells.

    Example
    -------

    Here is the example of annotation of dendritic cells based on their gene
    expressions. For annotation, we use data by Villani et al. (2017)[1] and
    marker genes by Franzén et al. (2019)[2].

    [1]  Villani, A. C., Satija, ... Jardine, L. (2017). Single-cell
         RNA-seq reveals new types of human blood dendritic cells, monocytes,
         and progenitors. Science, 356(6335).

    [2]  Oscar Franzén, Li-Ming Gan, Johan L M Björkegren, PanglaoDB:
         a web server for exploration of mouse and human single-cell RNA
         sequencing data, Database, Volume 2019, 2019.

    >>> import pandas as pd
    >>> from scanpy import AnnData
    >>> from scanpy.external.tl import annotator
    >>> import urllib.request
    >>>
    >>> # download data in a temporary directory
    >>> file_name_data, _ = urllib.request.urlretrieve(
    ...     "https://github.com/biolab/cell-annotation/releases/download/"
    ...     "0.1.0/DC_expMatrix_DCnMono.csv.gz")
    >>> file_name_markers, _ = urllib.request.urlretrieve(
    ...     "https://github.com/biolab/cell-annotation/releases/download/"
    ...     "0.1.0/panglao_gene_markers_human.csv.gz")
    >>>
    >>> # read data with pandas
    >>> df = pd.read_csv(file_name_data, compression="gzip").iloc[:, :-2]
    >>> df_markers = pd.read_csv(file_name_markers, compression="gzip")
    >>>
    >>> # transform data to AnnData
    >>> anndata = AnnData(df.values, var=df.columns.values)
    >>>
    >>> # run annotation
    >>> scores = annotator(anndata, df_markers, normalize=True)
    """

    try:
        from pointannotator.annotate_samples import AnnotateSamples
    except ImportError:
        raise ImportError(
            'Please install point-annotator: \n\t' 'pip install point-annotator'
        )

    data_df = pd.DataFrame(adata.X, columns=adata.var.values.flatten())
    if num_genes is None:
        num_genes = data_df.shape[1]
        warnings.warn(
            "The number of\norganisms genes is not provided. It is "
            "currently\nset to the number of genes of the dataset.\n"
            "If you want to change it please set `num_genes` "
            "parameter."
        )

    annotations = AnnotateSamples.annotate_samples(
        data_df,
        markers,
        num_genes=num_genes,
        return_nonzero_annotations=return_nonzero_annotations,
        p_threshold=p_threshold,
        p_value_fun=p_value_fun,
        z_threshold=z_threshold,
        scoring="scoring_" + scoring,
        normalize=normalize,
    )

    return AnnData(annotations.values, var=annotations.columns.values)
