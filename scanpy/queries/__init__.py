import pandas as pd


def mitochondrial_genes(host, org):
    """Mitochondrial gene symbols for specific organism through BioMart.

    Parameters
    ----------
    host : {{'www.ensembl.org', ...}}
        A valid BioMart host URL.
    org : {{'hsapiens', 'mmusculus'}}
        Organism to query. Currently available are human ('hsapiens') and mouse
        ('mmusculus').

    Returns
    -------
    A `pd.Index` containing mitochondrial gene symbols.
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
