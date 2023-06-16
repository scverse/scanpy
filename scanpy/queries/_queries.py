import collections.abc as cabc
from functools import singledispatch
from types import MappingProxyType
from typing import Any, Union, Optional, Iterable, Dict, Mapping

import pandas as pd
from anndata import AnnData

from ..get import rank_genes_groups_df
from .._utils import _doc_params


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


@_doc_params(doc_org=_doc_org, doc_host=_doc_host, doc_use_cache=_doc_use_cache)
def simple_query(
    org: str,
    attrs: Union[Iterable[str], str],
    *,
    filters: Optional[Dict[str, Any]] = None,
    host: str = "www.ensembl.org",
    use_cache: bool = False,
) -> pd.DataFrame:
    """\
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
    elif isinstance(attrs, cabc.Iterable):
        attrs = list(attrs)
    else:
        raise TypeError(f"attrs must be of type list or str, was {type(attrs)}.")
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


@_doc_params(doc_org=_doc_org, doc_host=_doc_host, doc_use_cache=_doc_use_cache)
def biomart_annotations(
    org: str,
    attrs: Iterable[str],
    *,
    host: str = "www.ensembl.org",
    use_cache: bool = False,
) -> pd.DataFrame:
    """\
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
    Dataframe containing annotations.

    Examples
    --------
    Retrieve genes coordinates and chromosomes

    >>> import scanpy as sc
    >>> annot = sc.queries.biomart_annotations(
            "hsapiens",
            ["ensembl_gene_id", "start_position", "end_position", "chromosome_name"],
        ).set_index("ensembl_gene_id")
    >>> adata.var[annot.columns] = annot
    """
    return simple_query(org=org, attrs=attrs, host=host, use_cache=use_cache)


@_doc_params(doc_org=_doc_org, doc_host=_doc_host, doc_use_cache=_doc_use_cache)
def gene_coordinates(
    org: str,
    gene_name: str,
    *,
    gene_attr: str = "external_gene_name",
    chr_exclude: Iterable[str] = (),
    host: str = "www.ensembl.org",
    use_cache: bool = False,
) -> pd.DataFrame:
    """\
    Retrieve gene coordinates for specific organism through BioMart.

    Parameters
    ----------
    {doc_org}
    gene_name
        The gene symbol (e.g. "hgnc_symbol" for human) for which to retrieve
        coordinates.
    gene_attr
        The biomart attribute the gene symbol should show up for.
    chr_exclude
        A list of chromosomes to exclude from query.
    {doc_host}
    {doc_use_cache}

    Returns
    -------
    Dataframe containing gene coordinates for the specified gene symbol.

    Examples
    --------
    >>> import scanpy as sc
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


@_doc_params(doc_org=_doc_org, doc_host=_doc_host, doc_use_cache=_doc_use_cache)
def mitochondrial_genes(
    org: str,
    *,
    attrname: str = "external_gene_name",
    host: str = "www.ensembl.org",
    use_cache: bool = False,
    chromosome: str = "MT",
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
    chromosome
        Mitochrondrial chromosome name used in BioMart for organism.

    Returns
    -------
    Dataframe containing identifiers for mitochondrial genes.

    Examples
    --------
    >>> import scanpy as sc
    >>> mito_gene_names = sc.queries.mitochondrial_genes("hsapiens")
    >>> mito_ensembl_ids = sc.queries.mitochondrial_genes("hsapiens", attrname="ensembl_gene_id")
    >>> mito_gene_names_fly = sc.queries.mitochondrial_genes("dmelanogaster", chromosome="mitochondrion_genome")
    """
    return simple_query(
        org,
        attrs=[attrname],
        filters={"chromosome_name": [chromosome]},
        host=host,
        use_cache=use_cache,
    )


@singledispatch
@_doc_params(doc_org=_doc_org)
def enrich(
    container: Union[Iterable[str], Mapping[str, Iterable[str]]],
    *,
    org: str = "hsapiens",
    gprofiler_kwargs: Mapping[str, Any] = MappingProxyType({}),
) -> pd.DataFrame:
    """\
    Get enrichment for DE results.

    This is a thin convenience wrapper around the very useful gprofiler_.

    This method dispatches on the first argument, leading to the following two
    signatures::

        enrich(container, ...)
        enrich(adata: AnnData, group, key: str, ...)

    Where::

        enrich(adata, group, key, ...) = enrich(adata.uns[key]["names"][group], ...)

    .. _gprofiler: https://pypi.org/project/gprofiler-official/#description

    Parameters
    ----------
    container
        Contains list of genes you'd like to search. If container is a `dict` all
        enrichment queries are made at once.
    adata
        AnnData object whose group will be looked for.
    group
        The group whose genes should be used for enrichment.
    key
        Key in `uns` to find group under.
    {doc_org}
    gprofiler_kwargs
        Keyword arguments to pass to `GProfiler.profile`, see gprofiler_. Some
        useful options are `no_evidences=False` which reports gene intersections,
        `sources=['GO:BP']` which limits gene sets to only GO biological processes and
        `all_results=True` which returns all results including the non-significant ones.
    **kwargs
        All other keyword arguments are passed to `sc.get.rank_genes_groups_df`. E.g.
        pval_cutoff, log2fc_min.

    Returns
    -------
    Dataframe of enrichment results.

    Examples
    --------
    Using `sc.queries.enrich` on a list of genes:

    >>> import scanpy as sc
    >>> sc.queries.enrich(['KLF4', 'PAX5', 'SOX2', 'NANOG'], org="hsapiens")
    >>> sc.queries.enrich({{'set1':['KLF4', 'PAX5'], 'set2':['SOX2', 'NANOG']}}, org="hsapiens")

    Using `sc.queries.enrich` on an :class:`anndata.AnnData` object:

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
    gprofiler_kwargs = dict(gprofiler_kwargs)
    for k in ["organism"]:
        if gprofiler_kwargs.get(k) is not None:
            raise ValueError(
                f"Argument `{k}` should be passed directly through `enrich`, "
                "not through `gprofiler_kwargs`"
            )
    return gprofiler.profile(container, organism=org, **gprofiler_kwargs)


@enrich.register(AnnData)
def _enrich_anndata(
    adata: AnnData,
    group: str,
    *,
    org: Optional[str] = "hsapiens",
    key: str = "rank_genes_groups",
    pval_cutoff: float = 0.05,
    log2fc_min: Optional[float] = None,
    log2fc_max: Optional[float] = None,
    gene_symbols: Optional[str] = None,
    gprofiler_kwargs: Mapping[str, Any] = MappingProxyType({}),
) -> pd.DataFrame:
    de = rank_genes_groups_df(
        adata,
        group=group,
        key=key,
        pval_cutoff=pval_cutoff,
        log2fc_min=log2fc_min,
        log2fc_max=log2fc_max,
        gene_symbols=gene_symbols,
    )
    if gene_symbols is not None:
        gene_list = list(de[gene_symbols].dropna())
    else:
        gene_list = list(de["names"].dropna())
    return enrich(gene_list, org=org, gprofiler_kwargs=gprofiler_kwargs)
