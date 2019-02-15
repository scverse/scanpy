from typing import Union, Optional
from collections.abc import Iterable
from ..utils import doc_params
from functools import singledispatch
from anndata import AnnData


_doc_org = """\
org : `str`
    Organism to query. Must be an organism in ensembl biomart. "hsapiens",
    "mmusculus", "drerio", etc.\
"""

_doc_host = """\
host : `str`, optional (default: "www.ensembl.org")
    A valid BioMart host URL. Alternative values include archive urls (like
    "grch37.ensembl.org") or regional mirrors (like "useast.ensembl.org").\
"""

_doc_use_cache = """\
use_cache : `bool`, optional (default: False)
    Whether pybiomart should use a cache for requests. Will create a
    `.pybiomart.sqlite` file in current directory if used.\
"""


@doc_params(doc_org=_doc_org, doc_host=_doc_host, doc_use_cache=_doc_use_cache)
def simple_query(org: str, attrs: Union[list, str], filters: dict = None,
                 host: str = "www.ensembl.org", use_cache: bool = False):
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
    elif isinstance(attrs, Iterable):
        attrs = list(attrs)
    else:
        raise TypeError("attrs must be of type list or str, was {}.".format(type(attrs)))
    try:
        from pybiomart import Server
    except ImportError:
        raise ImportError("You need to install the `pybiomart` module.")
    server = Server(host, use_cache=use_cache)
    dataset = server.marts["ENSEMBL_MART_ENSEMBL"].datasets[
        "{}_gene_ensembl".format(org)
    ]
    res = dataset.query(
        attributes=attrs, filters=filters, use_attr_names=True
    )
    return res


@doc_params(doc_org=_doc_org, doc_host=_doc_host, doc_use_cache=_doc_use_cache)
def biomart_annotations(org, attrs, host="www.ensembl.org", use_cache=False):
    """
    Retrieve gene annotations from ensembl biomart.

    Parameters
    ----------
    {doc_org}
    attrs : `List[str]`
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
def gene_coordinates(org, gene_name, gene_attr="external_gene_name", chr_exclude=[],
                     host="www.ensembl.org", use_cache=False):
    """
    Retrieve gene coordinates for specific organism through BioMart.

    Parameters
    ----------
    {doc_org}
    gene_name :
        The gene symbol (e.g. 'hgnc_symbol' for human) for which to retrieve
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
    """
    res = simple_query(org=org,
                       attrs=["chromosome_name", "start_position", "end_position"],
                       filters={gene_attr: gene_name}, host=host, use_cache=use_cache)
    return res[~res["chromosome_name"].isin(chr_exclude)]


@doc_params(doc_org=_doc_org, doc_host=_doc_host, doc_use_cache=_doc_use_cache)
def mitochondrial_genes(org, attrname="external_gene_name", host="www.ensembl.org",
                        use_cache=False):
    """\
    Mitochondrial gene symbols for specific organism through BioMart.

    Parameters
    ----------
    {doc_org}
    attrname : `str`, optional (default: "external_gene_name")
        Biomart attribute field to return. Possible values include
        "external_gene_name", "ensembl_gene_id", "hgnc_symbol", "mgi_symbol",
        and "zfin_id_symbol".
    {doc_host}
    {doc_use_cache}

    Returns
    -------
    An `np.array` containing identifiers for mitochondrial genes.
    """
    return simple_query(org,
                        attrs=[attrname],
                        filters={"chromosome_name": ["MT"]},
                        host=host,
                        use_cache=use_cache)


@singledispatch
def enrich(container, org: Optional[str] = None):
    """Get enrichment for DE results."""
    from gprofiler import gprofiler
    return gprofiler(container, org)


@enrich.register
def _enrich_anndata(adata: AnnData, group, key: str = "rank_genes_groups", org: Optional[str] = None):
    return enrich(adata.uns[key]["names"][group], org=org)
