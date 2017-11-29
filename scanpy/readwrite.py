# Author: Alex Wolf (http://falexwolf.de)
"""Reading and Writing
"""

import os
import sys
import h5py
import numpy as np
import pandas as pd
import time
from anndata import AnnData

from . import settings as sett
from . import logging as logg

avail_exts = {'csv', 'xlsx', 'txt', 'h5', 'soft.gz', 'txt.gz', 'mtx', 'tab', 'data'}
""" Available file formats for reading data. """


# --------------------------------------------------------------------------------
# Reading and Writing data files and AnnData objects
# --------------------------------------------------------------------------------


def read(filename, sheet=None, ext=None, delimiter=None,
         first_column_names=False, backup_url=None, return_dict=False,
         cache=None):
    """Read file and return AnnData object.

    To speed up reading, consider passing `cache=True`, which creates an hdf5
    cache file.

    Parameters
    ----------
    filename : str, None
        If the filename has no file extension, it is interpreted as a key for
        generating a filename via `sc.settings.writedir + filename +
        sc.settings.file_format_data`.  This is the same behavior as in
        `sc.read(filename, ...)`.
    sheet : str, optional (default: None)
        Name of sheet/table in hdf5 or Excel file.
    cache : bool or None, optional (default: False)
        If False, read from source, if True, read from fast 'h5' cache.
        If settings.recompute != 'read', this defaults to True.
    ext : str, optional (default: automatically inferred from filename)
        Extension that indicates the file type.
    delimiter : str, optional (default: None)
        Delimiter that separates data within text file. If None, will split at
        arbitrary number of white spaces, which is different from enforcing
        splitting at any single white space ' '.
    first_column_names : bool, optional (default: False)
        Assume the first column stores samplenames. This is only necessary if
        the sample are floats or integers: strings in the first column are
        automatically assumed to be row names.
    backup_url : str, optional (default: None)
        Retrieve the file from an URL if not present on disk.
    return_dict : bool, optional (default: False)
        Return dictionary instead of AnnData object.

    Returns
    -------
    adata : :class:`~sc.api.AnnData` object, dict
        Returns dict if `return_dict` was `True`.

    If a dict, the dict contains

    X : array-like or sparse matrix
        Data array for further processing, columns correspond to genes,
        rows correspond to samples.
    row_names : np.ndarray, optional
        Array storing the names of rows (experimental labels of samples).
    col_names : np.ndarray, optional
        Array storing the names of columns (gene names).
    """
    filename = str(filename)  # allow passing pathlib.Path objects
    if is_filename(filename):
        data = read_file(filename, sheet, ext, delimiter,
                         first_column_names, backup_url, cache)
        if isinstance(data, dict):
            return data if return_dict else AnnData(data)
        elif isinstance(data, AnnData):
            return data._to_dict_fixed_width_arrays() if return_dict else data
        else:
            raise ValueError('Do not know how to process read data.')

    # generate filename and read to dict
    filekey = filename
    filename = sett.writedir + filekey + '.' + sett.file_format_data
    if not os.path.exists(filename):
        raise ValueError('Reading with filekey "{}" failed, the '
                         'inferred filename "{}" does not exist. '
                         'If you intended to provide a filename, either '
                         'use a filename ending on one of the available extensions {} '
                         'or pass the parameter `ext`.'
                         .format(filekey, filename, avail_exts))
    data = read_file_to_dict(filename, ext=sett.file_format_data)
    return data if return_dict else AnnData(data)


def read_10x_h5(filename, genome='mm10'):
    """Read hdf5 file with naming conventions of 10X Genomics.

    Parameters
    ----------
    filename : str
        Filename.
    genome : str, optional (default: 'mm10')
        Genome group in hdf5 file.

    Returns
    -------
    adata : :class:`~sc.api.AnnData`
        Annotated data matrix, where samples/cells are named by their barcode
        and variables/genes by gene name. The data is stored in adata.X, cell
        names in adata.smp_names and gene names in adata.var_names.
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
                            {'smp_names': dsets['barcodes'].astype(str)},
                            {'var_names': dsets['gene_names'].astype(str),
                             'gene_ids': dsets['genes'].astype(str)})
            logg.info(t=True)
            return adata
        except tables.NoSuchNodeError:
            raise Exception('Genome %s does not exist in this file.' % genome)
        except KeyError:
            raise Exception('File is missing one or more required datasets.')


def write(filename, adata, ext=None, compression=None, compression_opts=None):
    """Write AnnData objects to file.

    Parameters
    ----------
    filename : `str`
        If the filename has no file extension, it is interpreted as a key for
        generating a filename via `sc.settings.writedir + filename +
        sc.settings.file_format_data`.  This is the same behavior as in
        `sc.read(filename, ...)`.
    adata : :class:`~scanpy.api.AnnData`
        Annotated data matrix.
    ext : {`None`, `'h5'`, `'csv'`, `'txt'`, `'npz'`} (default: `None`)
        File extension from wich to infer file format. If `None`, defaults to
        `sc.settings.file_format_data`.
    compression : None or {'gzip', 'lzf'}, optional (default: `None`)
        Only applies for hdf5 files. See
        http://docs.h5py.org/en/latest/high/dataset.html.
    compression_opts : int, optional (default: `None`)
        Only applies for hdf5 files. See
        See http://docs.h5py.org/en/latest/high/dataset.html.

    Examples
    --------
    Writing a usual numpy array.

    >>> adata = AnnData(
    >>>     data=np.array([[1, 0], [3, 0], [5, 6]]),
    >>>     smp={'row_names': ['name1', 'name2', 'name3'],
    >>>          'sanno1': ['cat1', 'cat2', 'cat2'],
    >>>          'sanno2': [2.1, 2.2, 2.3]},
    >>>     var={'vanno1': [3.1, 3.2]},
    >>>     uns={'sanno1_colors': ['#000000', '#FFFFFF'],
    >>>          'uns2': ['some annotation']})
    >>> assert pd.api.types.is_string_dtype(adata.smp['sanno1'])
    >>> write('./test.h5', adata)
    $ h5ls test.h5
    _data                    Dataset {3, 2}
    _smp                     Dataset {3}
    _var                     Dataset {2}
    sanno1_categories        Dataset {2}
    sanno1_colors            Dataset {2}
    uns2                     Dataset {1}
    $ h5ls -d test.h5
    ...
    _smp                     Dataset {3/3}
        Location:  1:1400
        Links:     1
        Storage:   42 logical bytes, 42 allocated bytes, 100.00% utilization
        Type:      struct {
                       "index"            +0    5-byte null-padded ASCII string
                       "sanno2"           +5    native double
                       "sanno1"           +13   native signed char
                   } 14 bytes
    ...

    Using compression

    >>> adata = AnnData(np.ones((10000, 1000)))
    >>> write('./test_compr.h5', adata, compression='gzip')
    >>> write('./test_no_compr.h5', adata)
    $ ls -lh test_no_compr.h5
    -rw-r--r--  1 alexwolf  staff    38M Nov 14 14:38 test_no_compr.h5
    $ ls -lh test_compr.h5
    -rw-r--r--  1 alexwolf  staff   161K Nov 14 15:15 test_compr.h5

    Writing a :class:`~scanpy.api.AnnData` that contains sparse data.

    >>> from scipy.sparse import csr_matrix
    >>> adata = AnnData(
    >>>     data=csr_matrix([[1, 0], [3, 0], [5, 6]]),
    >>>     smp={'row_names': ['name1', 'name2', 'name3'],
    >>>          'sanno1': ['cat1', 'cat2', 'cat2'],
    >>>          'sanno2': [2.1, 2.2, 2.3]},
    >>>     var={'vanno1': [3.1, 3.2]},
    >>>     uns={'sanno1_colors': ['#000000', '#FFFFFF'],
    >>>          'uns2_sparse': csr_matrix([[1, 0], [3, 0]])})
    $ h5ls test.h5
    _data_csr_data           Dataset {4}
    _data_csr_indices        Dataset {4}
    _data_csr_indptr         Dataset {4}
    _data_csr_shape          Dataset {2}
    _smp                     Dataset {3}
    _var                     Dataset {2}
    sanno1_categories        Dataset {2}
    sanno1_colors            Dataset {2}
    uns2_sparse_csr_data     Dataset {2}
    uns2_sparse_csr_indices  Dataset {2}
    uns2_sparse_csr_indptr   Dataset {3}
    uns2_sparse_csr_shape    Dataset {2}
    """
    filename = str(filename)  # allow passing pathlib.Path objects
    if is_filename(filename):
        filename = filename
        ext_ = is_filename(filename, return_ext=True)
        if ext is None:
            ext = ext_
        elif ext != ext_:
            raise ValueError('It suffices to provide the file type by '
                             'providing a proper extension to the filename.'
                             'One of "txt", "csv", "h5" or "npz".')
    else:
        key = filename
        ext = sett.file_format_data if ext is None else ext
        filename = get_filename_from_key(key, ext)
    write_anndata_to_file(filename, adata, ext=ext,
                          compression=compression,
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
    if not asheader:
        sett.m(verbosity, 'reading params file', filename)
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


def read_file(filename, sheet=None, ext=None, delimiter=None, first_column_names=None,
              backup_url=None, cache=None, suppress_cache_warning=False):
    """Read file and return data dictionary.

    To speed up reading and save storage space, this creates an hdf5 file if
    it's not present yet.

    Parameters
    ----------
    filename : str, Path
        Filename of data file.
    sheet : str, optional
        Name of sheet in Excel file.
    cache : bool or None (default: None)
        If false, read from source, if true, read from fast 'h5' cache.
    ext : str, optional (default: inferred from filename)
        Extension that indicates the file type.
    delimiter : str, optional
        Separator that separates data within text file. If None, will split at
        arbitrary number of white spaces, which is different from enforcing
        splitting at single white space ' '.
    first_column_names : bool, optional
        Assume the first column stores samplenames. Is unnecessary if the sample
        names are not floats or integers, but strings.
    backup_url : str
        URL for download of file in case it's not present.

    Returns
    -------
    ddata : dict containing
        X : np.ndarray
            Data array for further processing, columns correspond to genes,
            rows correspond to samples.
        row_names : np.ndarray
            Array storing the names of rows (experimental labels of samples).
        col_names : np.ndarray
            Array storing the names of columns (gene names).

    If sheet is unspecified and an h5 or xlsx file is read, the dict
    contains all sheets instead.
    """
    filename = str(filename)  # allow passing pathlib.Path objects
    if ext is not None and ext not in avail_exts:
        raise ValueError('Please provide one of the available extensions.\n'
                         + avail_exts)
    else:
        ext = is_filename(filename, return_ext=True)
    is_present = check_datafile_present_and_download(filename,
                                                     backup_url=backup_url)
    if not is_present: logg.m('... did not find original file', filename)
    # read hdf5 files
    if ext == 'h5':
        if sheet is None:
            return read_file_to_dict(filename, ext=sett.file_format_data)
        else:
            logg.msg('reading sheet', sheet, 'from file', filename, v=4)
            return _read_hdf5_single(filename, sheet)
    # read other file types
    filename_stripped = filename.lstrip('./')
    if filename_stripped.startswith('data/'):
        filename_stripped = filename_stripped[5:]
    fast_ext = sett.file_format_data if sett.file_format_data in {'h5', 'npz'} else 'h5'
    filename_cache = (sett.writedir + 'data/'
                      + filename_stripped.replace('.' + ext, '.' + fast_ext))
    cache = not sett.recompute == 'read' if cache is None else cache
    if cache and os.path.exists(filename_cache):
        ddata = read_file_to_dict(filename_cache, sett.file_format_data,
                                  cache_warning=True)
    else:
        if not is_present:
            raise FileNotFoundError(
                'Cannot read original data file {}, is not present.'
                .format(filename))
        logg.msg('reading', filename, v=4)
        if not cache and not suppress_cache_warning:
            logg.warn('This might be very slow. Consider passing `cache=True`, '
                      'which enables much faster reading from a cache file.')
        # do the actual reading
        if ext == 'xlsx' or ext == 'xls':
            if sheet is None:
                ddata = read_file_to_dict(filename, ext=ext)
            else:
                ddata = _read_excel(filename, sheet)
        elif ext == 'mtx':
            ddata = _read_mtx(filename)
        elif ext == 'csv':
            ddata = read_txt(filename, delimiter=',',
                             first_column_names=first_column_names)
        elif ext in {'txt', 'tab', 'data'}:
            if ext == 'data':
                logg.msg('... assuming ".data" means tab or white-space separated text file', v=3)
                logg.hint('--> change this by passing `ext` to sc.read')
            ddata = read_txt(filename, delimiter, first_column_names)
        elif ext == 'soft.gz':
            ddata = _read_softgz(filename)
        elif ext == 'txt.gz':
            sys.exit('TODO: implement similar to read_softgz')
        else:
            raise ValueError('Unkown extension {}.'.format(ext))
        if cache:
            logg.info('... writing an', sett.file_format_data,
                      'cache file to speedup reading next time')
            if not os.path.exists(os.path.dirname(filename_cache)):
                os.makedirs(os.path.dirname(filename_cache))
            # write for faster reading when calling the next time
            write_anndata_to_file(filename_cache, AnnData(ddata), sett.file_format_data)
    return ddata


def _read_mtx(filename, return_dict=True, dtype='float32'):
    """Read mtx file.
    """
    from scipy.io import mmread
    # could be rewritten accounting for dtype to be more performant
    X = mmread(filename).astype(dtype)
    from scipy.sparse import csr_matrix
    X = csr_matrix(X)
    logg.m('... did not find row_names or col_names')
    if return_dict:
        return {'X': X}
    else:
        return AnnData(X)


def read_txt(filename, delimiter=None, first_column_names=None, dtype='float32',
             return_dict=True):
    """Read text file and eturn data dictionary or AnnData object.

    Parameters
    ----------
    filename : str, Path
        Filename of data file.
    delimiter : str, optional
        Separator that separates data within text file. If None, will split at
        arbitrary number of white spaces, which is different from enforcing
        splitting at single white space ' '.
    first_column_names : bool, optional
        Assume the first column stores samplenames.
    return_dict : bool, optional (default: True)
        Return a dict, not an AnnData, if True.

    Returns
    -------
    data : An AnnData or a dict storing the data.

    If a dict, it contains
        X : np.ndarray
            Data array for further processing, columns correspond to genes,
            rows correspond to samples.
        row_names : np.ndarray
            Array storing the names of rows (experimental labels of samples).
        col_names : np.ndarray
            Array storing the names of columns (gene names).
    """
    filename = str(filename)  # allow passing pathlib.Path objects
    header = ''
    data = []
    length = -1
    f = open(filename)
    col_names = []
    row_names = []
    # read header and column names
    for line in f:
        if line.startswith('#'):
            header += line
        else:
            if delimiter is not None and delimiter not in line:
                logg.warn('did not find delimiter "{}" in first line'.format(delimiter))
            line_list = line.split(delimiter)
            if not is_float(line_list[0]):
                col_names = line_list
                logg.msg('    assuming first line in file stores column names', v=4)
            else:
                if not is_float(line_list[0]) or first_column_names:
                    first_column_names = True
                    row_names.append(line_list[0])
                    data.append(np.array(line_list[1:], dtype=dtype))
                else:
                    data.append(np.array(line_list, dtype=dtype))
            break
    if not col_names:
        # try reading col_names from the last comment line
        if len(header) > 0:
            logg.msg('    assuming last comment line stores variable names', v=4)
            col_names = np.array(header.split('\n')[-2].strip('#').split())
        # just numbers as col_names
        else:
            logg.msg('    did not find column names in file', v=4)
            col_names = np.arange(len(data[0])).astype(str)
    col_names = np.array(col_names, dtype=str)
    # read another line to check if first column contains row names or not
    if first_column_names is None: first_column_names = False
    for line in f:
        line_list = line.split(delimiter)
        if first_column_names or not is_float(line_list[0]):
            logg.msg('    assuming first column in file stores row names', v=4)
            first_column_names = True
            row_names.append(line_list[0])
            data.append(np.array(line_list[1:], dtype=dtype))
        else:
            data.append(np.array(line_list, dtype=dtype))
        break
    # if row names are just integers
    if len(data) > 1 and data[0].size != data[1].size:
        logg.msg('    assuming first row stores column names and first column row names', v=4)
        first_column_names = True
        col_names = np.array(data[0]).astype(int).astype(str)
        row_names.append(data[1][0].astype(int).astype(str))
        data = [data[1][1:]]
    # parse the file
    for line in f:
        line_list = line.split(delimiter)
        if first_column_names:
            row_names.append(line_list[0])
            data.append(np.array(line_list[1:], dtype=dtype))
        else:
            data.append(np.array(line_list, dtype=dtype))
    logg.msg('    read data into list of lists', t=True, v=4)
    # transfrom to array, this takes a long time and a lot of memory
    # but it's actually the same thing as np.genfromtext does
    # - we don't use the latter as it would involve another slicing step
    #   in the end, to separate row_names from float data, slicing takes
    #   a lot of memory and cpu time
    if data[0].size != data[-1].size:
        raise ValueError('length of first line {} is different from length of last line {}'
                         .format(data[0].size, data[-1].size))
    data = np.array(data, dtype=dtype)
    logg.msg('    constructed array from list of list', t=True, v=4)
    # transform row_names
    if not row_names:
        row_names = np.arange(len(data)).astype(str)
        logg.msg('    did not find row names in file', v=4)
    else:
        row_names = np.array(row_names)
        for iname, name in enumerate(row_names):
            row_names[iname] = name.strip('"')
    # adapt col_names if necessary
    if col_names.size > data.shape[1]:
        col_names = col_names[1:]
    for iname, name in enumerate(col_names):
        col_names[iname] = name.strip('"')
    if return_dict:
        return {'X': data, 'col_names': col_names, 'row_names': row_names}
    else:
        return AnnData(data, smp={'smp_names': row_names}, var={'var_names': col_names})


def _read_hdf5_single(filename, key=''):
    """Read a single dataset from an hdf5 file.

    See also function read_file_to_dict(), which reads all keys of the hdf5 file.

    Parameters
    ----------
    filename : str, Path
        Filename of data file.
    key : str, optional
        Name of dataset in the file. If not specified, shows available keys but
        raises an Error.

    Returns
    -------
    ddata : dict containing
        X : np.ndarray
            Data array for further processing, columns correspond to genes,
            rows correspond to samples.
        row_names : np.ndarray
            Array storing the names of rows (experimental labels of samples).
        col_names : np.ndarray
            Array storing the names of columns (gene names).
    """
    filename = str(filename)  # allow passing pathlib.Path objects
    with h5py.File(filename, 'r') as f:
        # the following is necessary in Python 3, because only
        # a view and not a list is returned
        keys = [k for k in f.keys()]
        if key == '':
            raise ValueError('The file ' + filename +
                             ' stores the following sheets:\n' + str(keys) +
                             '\n Call read/read_hdf5 with one of them.')
        # fill array
        X = f[key][()]
        if X.dtype.kind == 'S':
            X = X.astype(str)
        # init dict
        ddata = {'X': X}
        # try to find row and column names
        for iname, name in enumerate(['row_names', 'col_names']):
            if name in keys:
                ddata[name] = f[name][()]
            elif key + '_' + name.replace('_', '') in keys:
                ddata[name] = f[key + '_' + name.replace('_', '')][()]
            else:
                ddata[name] = np.arange(X.shape[0 if name == 'row_names' else 1])
                if key == 'X':
                    logg.info('did not find', name, 'in', filename)
            ddata[name] = ddata[name].astype(str)
            if X.ndim == 1:
                break
    return ddata


def _read_excel(filename, sheet):
    """Read excel file and return data dictionary.

    Parameters
    ----------
    filename : str, Path
        Filename to read from.
    sheet : str
        Name of sheet in Excel file.

    Returns
    -------
    ddata : dict containing
        X : np.ndarray
            Data array for further processing, columns correspond to genes,
            rows correspond to samples.
        row_names : np.ndarray
            Array storing the names of rows (experimental labels of samples).
        col_names : np.ndarray
            Array storing the names of columns (gene names).
    """
    filename = str(filename)  # allow passing pathlib.Path objects
    # rely on pandas for reading an excel file
    try:
        from pandas import read_excel
        df = read_excel(filename, sheet)
    except Exception as e:
        # in case this raises an error using Python 2.7
        print('if on Python 2.7 '
              'try installing xlrd using "sudo pip install xlrd"')
        raise e
    return ddata_from_df(df)


def _read_softgz(filename):
    """Read a SOFT format data file.

    The SOFT format is documented here
    http://www.ncbi.nlm.nih.gov/geo/info/soft2.html.

    Returns
    -------
    ddata : dict, containing
        X : np.ndarray
            A d × n array of gene expression values.
        col_names : np.ndarray
            A list of gene identifiers of length d.
        row_names : np.ndarray
            A list of sample identifiers of length n.
        groups : np.ndarray
            A list of sample desriptions of length n.

    Note
    ----
    The function is based on a script by Kerby Shedden.
    http://dept.stat.lsa.umich.edu/~kshedden/Python-Workshop/gene_expression_comparison.html
    """
    filename = str(filename)  # allow passing pathlib.Path objects
    import gzip
    with gzip.open(filename) as file:
        # The header part of the file contains information about the
        # samples. Read that information first.
        samples_info = {}
        for line in file:
            line = line.decode("utf-8")
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
        sample_names = file.readline().decode("utf-8").split("\t")
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
            line = line.decode("utf-8")
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
    smp = np.zeros((len(row_names),), dtype=[('smp_names', 'S21'), ('groups', 'S21')])
    smp['smp_names'] = sample_names
    smp['groups'] = groups
    var = np.zeros((len(gene_names),), dtype=[('var_names', 'S21')])
    var['var_names'] = gene_names
    ddata = {'X': X, 'smp': smp, 'var': var}
    return ddata


# -------------------------------------------------------------------------------
# Reading and writing for dictionaries
# -------------------------------------------------------------------------------


def read_file_to_dict(filename, ext='h5', cache_warning=False):
    """Read file and return dict with keys.

    The recommended format for this is hdf5.

    If reading from an Excel file, key names correspond to sheet names.

    Parameters
    ----------
    filename : str, Path
        Filename of data file.
    ext : {'h5', 'xlsx'}, optional
        Choose file format. Excel is much slower.

    Returns
    -------
    d : dict
    """
    filename = str(filename)  # allow passing pathlib.Path objects
    if cache_warning: logg.info('reading cache', filename)
    else: logg.msg('reading', filename, v=4)
    d = {}
    if ext == 'h5':
        with h5py.File(filename, 'r') as f:
            for key in f.keys():
                # the '()' means 'read everything' (by contrast, ':' only works
                # if not reading a scalar type)
                value = f[key][()]
                key, value = postprocess_reading(key, value)
                d[key] = value
    elif ext == 'npz':
        d_read = np.load(filename)
        for key, value in d_read.items():
            key, value = postprocess_reading(key, value)
            d[key] = value
    else:
        raise ValueError('`ext` need to be "h5" or "npz"')
    csr_keys = [key.replace('_csr_data', '')
                 for key in d if '_csr_data' in key]
    for key in csr_keys: d = load_sparse_csr(d, key=key)
    return d


def postprocess_reading(key, value):
    if value.ndim == 1 and len(value) == 1:
        value = value[0]
    if value.dtype.kind == 'S':
        value = value.astype(str)
        # recover a dictionary that has been stored as a string
        if len(value) > 0:
            if value[0] == '{' and value[-1] == '}': value = eval(value)
    if (key != 'smp' and key != 'var' and key != '_smp' and key != '_var'
        and not isinstance(value, dict) and value.dtype.names is not None):
        # TODO: come up with a better way of solving this, see also below
        new_dtype = [((dt[0], 'U{}'.format(int(int(dt[1][2:])/4)))
                      if dt[1][1] == 'S' else dt) for dt in value.dtype.descr]
        value = value.astype(new_dtype)
    return key, value


def preprocess_writing(key, value):
    if isinstance(value, dict):
        # hack for storing dicts
        value = np.array([str(value)])
    else:
        value = np.array(value)
        if value.ndim == 0:
            value = np.array([value])
    # some output about the data to write
    logg.m(key, type(value),
           value.dtype, value.dtype.kind, value.shape,
           v=6)
    # make sure string format is chosen correctly
    if value.dtype.kind == 'U':
        value = value.astype(np.string_)
    return key, value


def write_anndata_to_file(filename, adata, ext='h5',
                          compression=None, compression_opts=None):
    """Write dictionary to file.

    Values need to be np.arrays or transformable to numpy arrays.

    Parameters
    ----------
    filename : str, Path
        Filename of data file.
    d : dict
        Dictionary storing keys with np.ndarray-like data or scalars.
    ext : string
        Determines file type, allowed are 'h5' (hdf5),
        'xlsx' (Excel) [or 'csv' (comma separated value file)].
    """
    filename = str(filename)  # allow passing pathlib.Path objects
    directory = os.path.dirname(filename)
    if not os.path.exists(directory):
        logg.info('creating directory', directory + '/', 'for saving output files')
        os.makedirs(directory)
    # output the following at warning level, it's very important for the users
    if ext in {'h5', 'npz'}:
        logg.msg('writing', filename, v=4)
        d = adata._to_dict_fixed_width_arrays()
    else:
        d = adata._to_dict_dataframes()
    from scipy.sparse import issparse
    d_write = {}
    for key, value in d.items():
        if issparse(value):
            for k, v in save_sparse_csr(value, key=key).items():
                d_write[k] = v
        elif not isinstance(value, pd.DataFrame):
            key, value = preprocess_writing(key, value)
            d_write[key] = value
        else:
            d_write[key] = value
    # now open the file
    if ext == 'h5':
        if issparse(value):
            for k, v in save_sparse_csr(value, key=key).items():
                d_write[k] = v
        with h5py.File(filename, 'w') as f:
            for key, value in d_write.items():
                try:
                    # ignore arrays with empty dtypes
                    if value.dtype.descr:
                        f.create_dataset(key, data=value,
                                         compression=compression,
                                         compression_opts=compression_opts)
                except TypeError:
                    # try writing it as byte strings
                    try:
                        if value.dtype.names is None:
                            f.create_dataset(key, data=value.astype('S'),
                                             compression=compression,
                                             compression_opts=compression_opts)
                        else:
                            new_dtype = [(dt[0], 'S{}'.format(int(dt[1][2:])*4))
                                         for dt in value.dtype.descr]
                            f.create_dataset(key, data=value.astype(new_dtype),
                                             compression=compression,
                                             compression_opts=compression_opts)
                    except Exception as e:
                        logg.info(str(e))
                        logg.warn('Could not save field with key = "{}" to h5 file.'
                                  .format(key))
    elif ext == 'npz':
        np.savez(filename, **d_write)
    elif ext == 'csv' or ext == 'txt':
        # here this is actually a directory that corresponds to the
        # single hdf5 file
        dirname = filename.replace('.' + ext, '/')
        # write the following at warning level, it's very important for the users
        logg.info('writing', ext, 'files to', dirname)
        if not os.path.exists(dirname): os.makedirs(dirname)
        if not os.path.exists(dirname + 'uns'): os.makedirs(dirname + 'uns')
        not_yet_raised_data_graph_warning = True
        for key, value in d_write.items():
            if key.startswith('data_graph'):
                if not_yet_raised_data_graph_warning:
                    logg.warn('Omitting to write neighborhood graph (`adata.uns[\'data_graph...\']`).')
                    not_yet_raised_data_graph_warning = False
                continue
            filename = dirname
            if key not in {'data', 'var', 'smp', 'smpm', 'varm'}:
                filename += 'uns/'
            filename += key + '.' + ext
            df = value
            if not isinstance(value, pd.DataFrame):
                if np.ndim(value) == 0: value = value[None]
                try:
                    df = pd.DataFrame(value)
                except:
                    logg.warn('Omitting to write \'{}\'.'.format(key))
                    continue
            df.to_csv(filename, sep=(' ' if ext == 'txt' else ','),
                      header=True if key != 'data' else False,
                      index=True if key in {'smp', 'var'} else False)


# -------------------------------------------------------------------------------
# Type conversion
# -------------------------------------------------------------------------------


def save_sparse_csr(X, key='X'):
    from scipy.sparse.csr import csr_matrix
    X = csr_matrix(X)
    key_csr = key + '_csr'
    return {key_csr + '_data': X.data,
            key_csr + '_indices': X.indices,
            key_csr + '_indptr': X.indptr,
            key_csr + '_shape': np.array(X.shape)}


def load_sparse_csr(d, key='X'):
    from scipy.sparse.csr import csr_matrix
    key_csr = key + '_csr'
    d[key] = csr_matrix((d[key_csr + '_data'],
                         d[key_csr + '_indices'],
                         d[key_csr + '_indptr']),
                        shape=d[key_csr + '_shape'])
    del d[key_csr + '_data']
    del d[key_csr + '_indices']
    del d[key_csr + '_indptr']
    del d[key_csr + '_shape']
    return d


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
    ext = sett.file_format_data if ext is None else ext
    filename = sett.writedir + key + '.' + ext
    return filename


def ddata_from_df(df):
    """Write pandas.dataframe to ddata dictionary.
    """
    ddata = {
        'X': df.values[:, 1:].astype(float),
        'row_names': df.iloc[:, 0].values.astype(str),
        'col_names': np.array(df.columns[1:], dtype=str)}
    return ddata


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


def is_filename(filename, return_ext=False):
    """Check whether the argument is a filename."""
    for ext in avail_exts:
        if filename.endswith('.' + ext):
            return ext if return_ext else True
    if return_ext:
        raise ValueError('"{}" does not end on a valid extension.\n'
                         'Please, provide one of the available extensions.\n{}'
                         .format(filename, avail_exts))
    else:
        return False


# -------------------------------------------------------------------------------
# Tests
# -------------------------------------------------------------------------------


# need data files to run this test
# def test_profile_memory():
#     print()
#     logg.print_memory_usage()
#     # filename = 'data/paul15/paul15.h5'
#     # url = 'http://falexwolf.de/data/paul15.h5'
#     # adata = read(filename, 'data.debatched', backup_url=url)
#     filename = 'data/maehr17/blood_counts_raw.data'
#     adata = read(filename, reread=True)
#     logg.print_memory_usage()
