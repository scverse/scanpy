"""Reading and Writing
"""

import os
import sys
import numpy as np
import time
from pathlib import Path
from anndata import AnnData, read_loom, \
    read_csv, read_excel, read_text, read_hdf, read_mtx
from anndata import read as read_h5ad

from . import settings
from . import logging as logg

# .gz and .bz2 suffixes are also allowed for text formats
text_exts = {'csv',
             'tsv', 'tab', 'data', 'txt'} # these four are all equivalent
avail_exts = {'anndata', 'xlsx',
              'h5', 'h5ad',
              'soft.gz', 'mtx'} | text_exts
"""Available file formats for reading data. """


# --------------------------------------------------------------------------------
# Reading and Writing data files and AnnData objects
# --------------------------------------------------------------------------------


def read(filename, backed=False, sheet=None, ext=None, delimiter=None,
         first_column_names=False, backup_url=None, cache=False):
    """Read file and return :class:`~scanpy.api.AnnData` object.

    To speed up reading, consider passing `cache=True`, which creates an hdf5
    cache file.

    Parameters
    ----------
    filename : `str`
        If the filename has no file extension, it is interpreted as a key for
        generating a filename via `sc.settings.writedir + filename +
        sc.settings.file_format_data`.  This is the same behavior as in
        `sc.read(filename, ...)`.
    backed : {`False`, `True`, 'r', 'r+'}, optional (default: `False`)
        Load :class:`~scanpy.api.AnnData` in `backed` mode instead of fully
        loading it into memory (`memory` mode). Only applies to `.h5ad` files.
        `True` and 'r' are equivalent. If you want to modify backed attributes
        of the AnnData object, you need to choose 'r+'.
    sheet : `str`, optional (default: `None`)
        Name of sheet/table in hdf5 or Excel file.
    cache : `bool`, optional (default: `False`)
        If `False`, read from source, if `True`, read from fast 'h5ad' cache.
    ext : `str`, optional (default: `None`)
        Extension that indicates the file type. If `None`, uses extension of
        filename.
    delimiter : `str`, optional (default: `None`)
        Delimiter that separates data within text file. If `None`, will split at
        arbitrary number of white spaces, which is different from enforcing
        splitting at any single white space ' '.
    first_column_names : `bool`, optional (default: `False`)
        Assume the first column stores row names. This is only necessary if
        these are not strings: strings in the first column are automatically
        assumed to be row names.
    backup_url : `str`, optional (default: `None`)
        Retrieve the file from an URL if not present on disk.

    Returns
    -------
    adata : :class:`~scanpy.api.AnnData`
    """
    filename = str(filename)  # allow passing pathlib.Path objects
    if is_valid_filename(filename):
        return _read(filename, backed=backed, sheet=sheet, ext=ext,
                     delimiter=delimiter, first_column_names=first_column_names,
                     backup_url=backup_url, cache=cache)
    # generate filename and read to dict
    filekey = filename
    filename = settings.writedir + filekey + '.' + settings.file_format_data
    if not os.path.exists(filename):
        raise ValueError('Reading with filekey "{}" failed, the '
                         'inferred filename "{}" does not exist. '
                         'If you intended to provide a filename, either '
                         'use a filename ending on one of the available extensions {} '
                         'or pass the parameter `ext`.'
                         .format(filekey, filename, avail_exts))
    return read_h5ad(filename, backed=backed)


def read_10x_h5(filename, genome='mm10'):
    """Read 10x-Genomics-formatted hdf5 file.

    Parameters
    ----------
    filename : `str`
        Filename.
    genome : `str`, optional (default: 'mm10')
        Genome group in hdf5 file.

    Returns
    -------
    adata : :class:`~scanpy.api.AnnData`
        Annotated data matrix, where obsevations/cells are named by their
        barcode and variables/genes by gene name. The data matrix is stored in
        `adata.X`, cell names in `adata.obs_names` and gene names in
        `adata.var_names`. The gene IDs are stored in `adata.obs['gene_ids']`.
    """
    logg.info('reading', filename, r=True, end=' ')
    import tables
    with tables.open_file(filename, 'r') as f:
        try:
            dsets = {}
            for node in f.walk_nodes('/' + genome, 'Array'):
                dsets[node.name] = node.read()
            # AnnData works with csr matrices
            # 10x stores the transposed data, so we do the transposition right away
            from scipy.sparse import csr_matrix
            M, N = dsets['shape']
            data = dsets['data']
            if dsets['data'].dtype == np.dtype('int32'):
                data = dsets['data'].view('float32')
                data[:] = dsets['data']
            matrix = csr_matrix((data, dsets['indices'], dsets['indptr']),
                                shape=(N, M))
            # the csc matrix is automatically the transposed csr matrix
            # as scanpy expects it, so, no need for a further transpostion
            adata = AnnData(matrix,
                            {'obs_names': dsets['barcodes'].astype(str)},
                            {'var_names': dsets['gene_names'].astype(str),
                             'gene_ids': dsets['genes'].astype(str)})
            logg.info(t=True)
            return adata
        except tables.NoSuchNodeError:
            raise Exception('Genome %s does not exist in this file.' % genome)
        except KeyError:
            raise Exception('File is missing one or more required datasets.')


def write(filename, adata, ext=None, compression='gzip', compression_opts=None):
    """Write :class:`~scanpy.api.AnnData` objects to file.

    Parameters
    ----------
    filename : `str`
        If the filename has no file extension, it is interpreted as a key for
        generating a filename via `sc.settings.writedir + filename +
        sc.settings.file_format_data`.  This is the same behavior as in
        :func:`~scanpy.api.read`.
    adata : :class:`~scanpy.api.AnnData`
        Annotated data matrix.
    ext : {`None`, `'h5'`, `'csv'`, `'txt'`, `'npz'`} (default: `None`)
        File extension from wich to infer file format. If `None`, defaults to
        `sc.settings.file_format_data`.
    compression : {`None`, 'gzip', 'lzf'}, optional (default: `'gzip'`)
        See http://docs.h5py.org/en/latest/high/dataset.html.
    compression_opts : `int`, optional (default: `None`)
        See http://docs.h5py.org/en/latest/high/dataset.html.
    """
    filename = str(filename)  # allow passing pathlib.Path objects
    if is_valid_filename(filename):
        filename = filename
        ext_ = is_valid_filename(filename, return_ext=True)
        if ext is None:
            ext = ext_
        elif ext != ext_:
            raise ValueError('It suffices to provide the file type by '
                             'providing a proper extension to the filename.'
                             'One of "txt", "csv", "h5" or "npz".')
    else:
        key = filename
        ext = settings.file_format_data if ext is None else ext
        filename = get_filename_from_key(key, ext)
    if ext == 'csv':
        adata.write_csvs(filename)
    else:
        adata.write(filename, compression=compression,
                    compression_opts=compression_opts)


# -------------------------------------------------------------------------------
# Reading and writing parameter files
# -------------------------------------------------------------------------------


def read_params(filename, asheader=False, verbosity=0):
    """Read parameter dictionary from text file.

    Assumes that parameters are specified in the format:
        par1 = value1
        par2 = value2

    Comments that start with '#' are allowed.

    Parameters
    ----------
    filename : str, Path
        Filename of data file.
    asheader : bool, optional
        Read the dictionary from the header (comment section) of a file.

    Returns
    -------
    params : dict
        Dictionary that stores parameters.
    """
    filename = str(filename)  # allow passing pathlib.Path objects
    from collections import OrderedDict
    params = OrderedDict([])
    for line in open(filename):
        if '=' in line:
            if not asheader or line.startswith('#'):
                line = line[1:] if line.startswith('#') else line
                key, val = line.split('=')
                key = key.strip()
                val = val.strip()
                params[key] = convert_string(val)
    return params


def write_params(filename, *args, **dicts):
    """Write parameters to file, so that it's readable by read_params.

    Uses INI file format.
    """
    if not os.path.exists(os.path.dirname(filename)):
        os.makedirs(os.path.dirname(filename))
    if len(args) == 1:
        d = args[0]
        with open(filename, 'w') as f:
            for key in d:
                f.write(key + ' = ' + str(d[key]) + '\n')
    else:
        with open(filename, 'w') as f:
            for k, d in dicts.items():
                f.write('[' + k + ']\n')
                for key, val in d.items():
                    f.write(key + ' = ' + str(val) + '\n')


def get_params_from_list(params_list):
    """Transform params list to dictionary.
    """
    params = {}
    for i in range(0, len(params_list)):
        if '=' not in params_list[i]:
            try:
                if not isinstance(params[key], list): params[key] = [params[key]]
                params[key] += [params_list[i]]
            except KeyError:
                raise ValueError('Pass parameters like `key1=a key2=b c d key3=...`.')
        else:
            key_val = params_list[i].split('=')
            key, val = key_val
            params[key] = convert_string(val)
    return params


# -------------------------------------------------------------------------------
# Reading and Writing data files
# -------------------------------------------------------------------------------


def _read(filename, backed=False, sheet=None, ext=None, delimiter=None,
          first_column_names=None, backup_url=None, cache=False,
          suppress_cache_warning=False):
    if ext is not None and ext not in avail_exts:
        raise ValueError('Please provide one of the available extensions.\n'
                         + avail_exts)
    else:
        ext = is_valid_filename(filename, return_ext=True)
    is_present = check_datafile_present_and_download(filename,
                                                     backup_url=backup_url)
    if not is_present: logg.msg('... did not find original file', filename)
    # read hdf5 files
    if ext in {'h5', 'h5ad'}:
        if sheet is None:
            return read_h5ad(filename, backed=backed)
        else:
            logg.msg('reading sheet', sheet, 'from file', filename, v=4)
            return read_hdf(filename, sheet)
    # read other file types
    filename_cache = (settings.cachedir + filename.lstrip(
        './').replace('/', '-').replace('.' + ext, '.h5ad'))
    if cache and os.path.exists(filename_cache):
        logg.info('... reading from cache file', filename_cache)
        adata = read_h5ad(filename_cache, backed=False)
    else:
        if not is_present:
            raise FileNotFoundError('Did not find file {}.'.format(filename))
        logg.msg('reading', filename, v=4)
        if not cache and not suppress_cache_warning:
            logg.hint('This might be very slow. Consider passing `cache=True`, '
                      'which enables much faster reading from a cache file.')
        # do the actual reading
        if ext == 'xlsx' or ext == 'xls':
            if sheet is None:
                raise ValueError(
                    'Provide `sheet` parameter when reading \'.xlsx\' files.')
            else:
                adata = read_excel(filename, sheet)
        elif ext == 'mtx':
            adata = read_mtx(filename)
        elif ext == 'csv':
            adata = read_csv(filename, first_column_names=first_column_names)
        elif ext in {'txt', 'tab', 'data', 'tsv'}:
            if ext == 'data':
                logg.msg('... assuming \'.data\' means tab or white-space '
                         'separated text file', v=3)
                logg.hint('change this by passing `ext` to sc.read')
            adata = read_text(filename, delimiter, first_column_names)
        elif ext == 'soft.gz':
            adata = _read_softgz(filename)
        else:
            raise ValueError('Unkown extension {}.'.format(ext))
        if cache:
            logg.info('... writing an', settings.file_format_data,
                      'cache file to speedup reading next time')
            if not os.path.exists(os.path.dirname(filename_cache)):
                os.makedirs(os.path.dirname(filename_cache))
            # write for faster reading when calling the next time
            adata.write(filename_cache)
    return adata


def _read_softgz(filename):
    """Read a SOFT format data file.

    The SOFT format is documented here
    http://www.ncbi.nlm.nih.gov/geo/info/soft2.html.

    Returns
    -------
    adata

    Note
    ----
    The function is based on a script by Kerby Shedden.
    http://dept.stat.lsa.umich.edu/~kshedden/Python-Workshop/gene_expression_comparison.html
    """
    filename = str(filename)  # allow passing pathlib.Path objects
    import gzip
    with gzip.open(filename, mode='rt') as file:
        # The header part of the file contains information about the
        # samples. Read that information first.
        samples_info = {}
        for line in file:
            if line.startswith("!dataset_table_begin"):
                break
            elif line.startswith("!subset_description"):
                subset_description = line.split("=")[1].strip()
            elif line.startswith("!subset_sample_id"):
                subset_ids = line.split("=")[1].split(",")
                subset_ids = [x.strip() for x in subset_ids]
                for k in subset_ids:
                    samples_info[k] = subset_description
        # Next line is the column headers (sample id's)
        sample_names = file.readline().strip().split("\t")
        # The column indices that contain gene expression data
        I = [i for i, x in enumerate(sample_names) if x.startswith("GSM")]
        # Restrict the column headers to those that we keep
        sample_names = [sample_names[i] for i in I]
        # Get a list of sample labels
        groups = [samples_info[k] for k in sample_names]
        # Read the gene expression data as a list of lists, also get the gene
        # identifiers
        gene_names, X = [], []
        for line in file:
            # This is what signals the end of the gene expression data
            # section in the file
            if line.startswith("!dataset_table_end"):
                break
            V = line.split("\t")
            # Extract the values that correspond to gene expression measures
            # and convert the strings to numbers
            x = [float(V[i]) for i in I]
            X.append(x)
            gene_names.append(V[1])
    # Convert the Python list of lists to a Numpy array and transpose to match
    # the Scanpy convention of storing samples in rows and variables in colums.
    X = np.array(X).T
    row_names = sample_names
    col_names = gene_names
    obs = np.zeros((len(row_names),), dtype=[('obs_names', 'S21'), ('groups', 'S21')])
    obs['obs_names'] = sample_names
    obs['groups'] = groups
    var = np.zeros((len(gene_names),), dtype=[('var_names', 'S21')])
    var['var_names'] = gene_names
    ddata = {'X': X, 'obs': obs, 'var': var}
    return AnnData(ddata)


# -------------------------------------------------------------------------------
# Type conversion
# -------------------------------------------------------------------------------


def is_float(string):
    """Check whether string is float.

    See also
    --------
    http://stackoverflow.com/questions/736043/checking-if-a-string-can-be-converted-to-float-in-python
    """
    try:
        float(string)
        return True
    except ValueError:
        return False


def is_int(string):
    """Check whether string is integer.
    """
    try:
        int(string)
        return True
    except ValueError:
        return False


def convert_bool(string):
    """Check whether string is boolean.
    """
    if string == 'True':
        return True, True
    elif string == 'False':
        return True, False
    else:
        return False, False


def convert_string(string):
    """Convert string to int, float or bool.
    """
    if is_int(string):
        return int(string)
    elif is_float(string):
        return float(string)
    elif convert_bool(string)[0]:
        return convert_bool(string)[1]
    elif string == 'None':
        return None
    else:
        return string


# -------------------------------------------------------------------------------
# Helper functions for reading and writing
# -------------------------------------------------------------------------------


def get_used_files():
    """Get files used by processes with name scanpy."""
    import psutil
    loop_over_scanpy_processes = (proc for proc in psutil.process_iter()
                                  if proc.name() == 'scanpy')
    filenames = []
    for proc in loop_over_scanpy_processes:
        try:
            flist = proc.open_files()
            for nt in flist:
                filenames.append(nt.path)
        # This catches a race condition where a process ends
        # before we can examine its files
        except psutil.NoSuchProcess as err:
            pass
    return set(filenames)


def wait_until_file_unused(filename):
    while (filename in get_used_files()):
        time.sleep(1)


def get_filename_from_key(key, ext=None):
    ext = settings.file_format_data if ext is None else ext
    filename = settings.writedir + key + '.' + ext
    return filename


def download_progress(count, blockSize, totalSize):
    percent = int(count*blockSize*100/totalSize)
    sys.stdout.write('\r' + '... %d%%' % percent)
    sys.stdout.flush()


def check_datafile_present_and_download(filename, backup_url=None):
    """Check whether the file is present, otherwise download.
    """
    if os.path.exists(filename): return True
    if backup_url is None: return False
    logg.info('try downloading from url\n' + backup_url + '\n' +
              '... this may take a while but only happens once')
    d = os.path.dirname(filename)
    if not os.path.exists(d):
        logg.info('creating directory', d + '/', 'for saving data')
        os.makedirs(d)
    from urllib.request import urlretrieve
    urlretrieve(backup_url, filename, reporthook=download_progress)
    logg.info('')
    return True


def is_valid_filename(filename, return_ext=False):
    """Check whether the argument is a filename."""
    ext = Path(filename).suffixes

    # cases for gzipped/bzipped text files
    if len(ext) == 2 and ext[0][1:] in text_exts and ext[1][1:] in ('gz', 'bz2'):
        return ext[0][1:] if return_ext else True
    elif len(ext) == 1 and ext[0][1:] in avail_exts:
        return ext[0][1:] if return_ext else True
    elif ''.join(ext) == '.soft.gz':
        return 'soft.gz' if return_ext else True
    else:
        if return_ext:
            raise ValueError('"{}" does not end on a valid extension.\n'
                             'Please, provide one of the available extensions.\n{}\n'
                             'Text files with .gz and .bz2 extensions are also supported.'
                             .format(filename, avail_exts))
        else:
            return False
