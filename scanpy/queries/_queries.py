from copy import copy
from functools import singledispatch
from collections import abc
from typing import Any, Union, Optional, Iterable, Dict

from anndata import AnnData
import pandas as pd

from ..utils import doc_params


_doc_org = """\
org
    Organism to query. Must be an organism in ensembl biomart. "hsapiens",
    "mmusculus", "drerio", etc.\
"""

_doc_host = """\
host
    A valid BioMart host URL. Alternative values include archive urls (like
    "grch37.ensembl.org") or regional mirrors (like "useast.ensembl.org").\
"""

_doc_use_cache = """\
use_cache
    Whether pybiomart should use a cache for requests. Will create a
    `.pybiomart.sqlite` file in current directory if used.\
"""


@doc_params(doc_org=_doc_org, doc_host=_doc_host, doc_use_cache=_doc_use_cache)
def simple_query(
    org: str,
    attrs: Union[Iterable[str], str],
    *,
    filters: Optional[Dict[str, Any]] = None,
    host: str = "www.ensembl.org",
    use_cache: bool = False
) -> pd.DataFrame:
    """
    A simple interface to biomart.

    Params
    ------
    {doc_org}
    attrs
        What you want returned.
    filters
        What you want to pick out.
    {doc_host}
    {doc_use_cache}
    """
    if isinstance(attrs, str):
        attrs = [attrs]
    elif isinstance(attrs, abc.Iterable):
        attrs = list(attrs)
    else:
        raise TypeError(
            "attrs must be of type list or str, was {}.".format(type(attrs))
        )
    try:
        from pybiomart import Server
    except ImportError:
        raise ImportError(
            "This method requires the `pybiomart` module to be installed."
        )
    server = Server(host, use_cache=use_cache)
    dataset = server.marts["ENSEMBL_MART_ENSEMBL"].datasets[
        "{}_gene_ensembl".format(org)
    ]
    res = dataset.query(attributes=attrs, filters=filters, use_attr_names=True)
    return res


@doc_params(doc_org=_doc_org, doc_host=_doc_host, doc_use_cache=_doc_use_cache)
def biomart_annotations(
    org: str,
    attrs: Iterable[str],
    *,
    host: str = "www.ensembl.org",
    use_cache: bool = False
) -> pd.DataFrame:
    """
    Retrieve gene annotations from ensembl biomart.

    Parameters
    ----------
    {doc_org}
    attrs
        Attributes to query biomart for.
    {doc_host}
    {doc_use_cache}

    Returns
    -------
    A `pd.DataFrame` containing annotations.

    Examples
    --------
    Retrieve genes coordinates and chromosomes

    >>> annot = sc.queries.biomart_annotations(
            "hsapiens",
            ["ensembl_gene_id", "start_position", "end_position", "chromosome_name"],
        ).set_index("ensembl_gene_id")
    >>> adata.var[annot.columns] = annot
    """
    return simple_query(org=org, attrs=attrs, host=host, use_cache=use_cache)


@doc_params(doc_org=_doc_org, doc_host=_doc_host, doc_use_cache=_doc_use_cache)
def gene_coordinates(
    org: str,
    gene_name: str,
    *,
    gene_attr: str = "external_gene_name",
    chr_exclude: Iterable[str] = (),
    host: str = "www.ensembl.org",
    use_cache: bool = False
) -> pd.DataFrame:
    """
    Retrieve gene coordinates for specific organism through BioMart.

    Parameters
    ----------
    {doc_org}
    gene_name :
        The gene symbol (e.g. "hgnc_symbol" for human) for which to retrieve
        coordinates.
    gene_attr : `str`, optional (default: "external_gene_name")
        The biomart attribute the gene symbol should show up for.
    chr_exclude :
        A list of chromosomes to exclude from query.
    {doc_host}
    {doc_use_cache}

    Returns
    -------
    A `pd.DataFrame` containing gene coordinates for the specified gene symbol.

    Examples
    --------
    >>> sc.queries.gene_coordinates("hsapiens", "MT-TF")
    """
    res = simple_query(
        org=org,
        attrs=["chromosome_name", "start_position", "end_position"],
        filters={gene_attr: gene_name},
        host=host,
        use_cache=use_cache,
    )
    return res[~res["chromosome_name"].isin(chr_exclude)]


@doc_params(doc_org=_doc_org, doc_host=_doc_host, doc_use_cache=_doc_use_cache)
def mitochondrial_genes(
    org: str,
    *,
    attrname: str = "external_gene_name",
    host: str = "www.ensembl.org",
    use_cache: bool = False
) -> pd.DataFrame:
    """\
    Mitochondrial gene symbols for specific organism through BioMart.

    Parameters
    ----------
    {doc_org}
    attrname
        Biomart attribute field to return. Possible values include
        "external_gene_name", "ensembl_gene_id", "hgnc_symbol", "mgi_symbol",
        and "zfin_id_symbol".
    {doc_host}
    {doc_use_cache}

    Returns
    -------
    An `pd.DataFrame` containing identifiers for mitochondrial genes.

    Examples
    --------
    >>> mito_gene_names = sc.queries.mitochondrial_genes("hsapiens")
    >>> mito_ensembl_ids = sc.queries.mitochondrial_genes("hsapiens", attrname="ensembl_gene_id")
    """
    return simple_query(
        org,
        attrs=[attrname],
        filters={"chromosome_name": ["MT"]},
        host=host,
        use_cache=use_cache,
    )


@singledispatch
def enrich(
    container: Iterable[str],
    *,
    org: str = "hsapiens",
    gprofiler_kwargs: dict = {}
) -> Optional[pd.DataFrame]:
    """
    Get enrichment for DE results.

    This is a thin convenience wrapper around the very useful
    `gprofiler <https://github.com/vals/python-gprofiler>`_.

    This method dispatches on the first argument, leading to the following two
    signatures::

        enrich(container, ...)
        enrich(adata: AnnData, group, key: str, ...)

    Where::

        enrich(adata, group, key, ...) = enrich(adata.uns[key]["names"][group], ...)

    Parameters
    ----------
    container
        Contains genes you'd like to search.
    adata
        AnnData object whose group will be looked for.
    group
        The group whose genes should be used for enrichment.
    key
        Key in `uns` to find group under.
    org
        Organism to query. Must be an organism in ensembl biomart. "hsapiens",
        "mmusculus", "drerio", etc.
    gprofiler_kwargs
        Keyword arguments to pass to `gprofiler`.

    Returns
    -------
    `pd.DataFrame` or `None`:
        Returns enrichment results. If nothing was found to be enriched,
        returns `None`.

    Examples
    --------
    Using `sc.queries.enrich` on a list of genes:

    >>> sc.queries.enrich(['Klf4', 'Pax5', 'Sox2', 'Nanog'], org="hsapiens")

    Using `sc.queries.enrich on an `AnnData` object:

    >>> pbmcs = sc.datasets.pbmc68k_reduced()
    >>> sc.tl.rank_genes_groups(pbmcs, "bulk_labels")
    >>> sc.queries.enrich(pbmcs, "CD34+")
    """
    try:
        from gprofiler import GProfiler
    except ImportError:
        raise ImportError(
            "This method requires the `gprofiler-official` module to be installed."
        )
    gprofiler = GProfiler(user_agent="scanpy", return_dataframe=True)
    gprofiler_kwargs = copy(gprofiler_kwargs)
    for k in ["organism"]:
        if gprofiler_kwargs.get(k) is not None:
            raise ValueError(
                "Argument `{}` should be passed directly through `enrich`, not"
                " through `gprofiler_kwargs`".format(k)
            )
    return gprofiler.profile(
        list(container),
        organism=org,
        **gprofiler_kwargs
    )


@enrich.register
def _enrich_anndata(
    adata: AnnData,
    group: str,
    *,
    org: Optional[str] = "hsapiens",
    key: str = "rank_genes_groups",
    pval_cutoff: float = 0.05,
    logfc_cutoff: Optional[float] = None,
    gene_symbols: Optional[str] = None,
    gprofiler_kwargs: dict = {}
) -> Optional[pd.DataFrame]:
    de = rank_genes_groups_df(
        adata,
        group=group,
        key=key,
        pval_cutoff=pval_cutoff,
        logfc_cutoff=logfc_cutoff,
        gene_symbols=gene_symbols
    )
    if gene_symbols is not None:
        gene_list = list(de[gene_symbols])
    else:
        gene_list = list(de["names"])
    return enrich(gene_list, org=org, gprofiler_kwargs=gprofiler_kwargs)


####### Utilities

def rank_genes_groups_df(
    adata: AnnData,
    group: str,  # Can this be something else?
    key: str = "rank_genes_groups",
    pval_cutoff : Optional[float] = None,
    logfc_cutoff : Optional[float] = None,
    gene_symbols : Optional[str] = None
):
    d = pd.DataFrame()
    for k in ['scores', 'names', 'logfoldchanges', 'pvals', 'pvals_adj']:
        d[k] = adata.uns["rank_genes_groups"][k][group]
    if pval_cutoff is not None:
        d = d[d["pvals_adj"] < pval_cutoff]
    if logfc_cutoff is not None:
        d = d[d["logfoldchanges"].abs() > logfc_cutoff]
    if gene_symbols is not None:
        d = d.join(adata.var[gene_symbols], on="names")
    return d