import pandas as pd
from .. import logging as logg


def mitochondrial_genes(host, org) -> pd.Index:
    """Mitochondrial gene symbols for specific organism through BioMart.

    Parameters
    ----------
    host : {{'www.ensembl.org', ...}}
        A valid BioMart host URL.
    org : {{'hsapiens', 'mmusculus', 'drerio'}}
        Organism to query. Currently available are human ('hsapiens'), mouse
        ('mmusculus') and zebrafish ('drerio').

    Returns
    -------
    A :class:`pandas.Index` containing mitochondrial gene symbols.
    """
    try:
        from bioservices import biomart
    except ImportError:
        raise ImportError(
            'You need to install the `bioservices` module.')
    from io import StringIO
    s = biomart.BioMart(host=host)

    # building query
    s.new_query()
    if org == 'hsapiens':
        s.add_dataset_to_xml('hsapiens_gene_ensembl')
        s.add_attribute_to_xml('hgnc_symbol')
    elif org == 'mmusculus':
        s.add_dataset_to_xml('mmusculus_gene_ensembl')
        s.add_attribute_to_xml('mgi_symbol')
    elif org == 'drerio':
        s.add_dataset_to_xml('drerio_gene_ensembl')
        s.add_attribute_to_xml('zfin_id_symbol')
    else:
        logg.msg('organism ', str(org), ' is unavailable', v=4, no_indent=True)
        return None
    s.add_attribute_to_xml('chromosome_name')
    xml = s.get_xml()

    # parsing mitochondrial gene symbols
    res = pd.read_csv(StringIO(s.query(xml)), sep='\t', header=None)
    res.columns = ['symbol', 'chromosome_name']
    res = res.dropna()
    res = res[res['chromosome_name'] == 'MT']
    res = res.set_index('symbol')
    res = res[~res.index.duplicated(keep='first')]

    return res.index


def gene_coordinates(host, org, gene, chr_exclude=[]) -> pd.DataFrame:
    """Retrieve gene coordinates for specific organism through BioMart.
    Parameters
    ----------
    host : {{'www.ensembl.org', ...}}
        A valid BioMart host URL. Can be used to control genome build.
    org : {{'hsapiens', 'mmusculus', 'drerio'}}
        Organism to query. Currently available are human ('hsapiens'), mouse
        ('mmusculus') and zebrafish ('drerio').
    gene :
        The gene symbol (e.g. 'hgnc_symbol' for human) for which to retrieve
        coordinates.
    chr_exclude :
        A list of chromosomes to exclude from query.
    Returns
    -------
    A `pd.DataFrame` containing gene coordinates for the specified gene symbol.
    """
    try:
        from bioservices import biomart
    except ImportError:
        raise ImportError(
            'You need to install the `bioservices` module.')
    from io import StringIO
    s = biomart.BioMart(host=host)

    # building query
    s.new_query()
    if org == 'hsapiens':
        s.add_dataset_to_xml('hsapiens_gene_ensembl')
        s.add_attribute_to_xml('hgnc_symbol')
    elif org == 'mmusculus':
        s.add_dataset_to_xml('mmusculus_gene_ensembl')
        s.add_attribute_to_xml('mgi_symbol')
    elif org == 'drerio':
        s.add_dataset_to_xml('drerio_gene_ensembl')
        s.add_attribute_to_xml('zfin_id_symbol')
    else:
        logg.msg('organism ', str(org), ' is unavailable', v=4, no_indent=True)
        return None
    s.add_attribute_to_xml('chromosome_name')
    s.add_attribute_to_xml('start_position')
    s.add_attribute_to_xml('end_position')
    xml = s.get_xml()

    # parsing gene coordinates
    res = pd.read_csv(StringIO(s.query(xml)), sep='\t', header=None)
    res.columns = ['symbol', 'chromosome_name', 'start', 'end']
    res = res.dropna()
    res = res[~res['chromosome_name'].isin(chr_exclude)]
    res = res.set_index('symbol')

    return res.loc[[gene], :]
