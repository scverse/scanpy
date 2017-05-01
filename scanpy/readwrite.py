# Author: F. Alex Wolf (http://falexwolf.de)
"""
Reading and Writing
"""

import os
import sys
import h5py
import numpy as np

from . import settings as sett

avail_exts = ['csv', 'xlsx', 'txt', 'h5', 'soft.gz', 'txt.gz', 'mtx', 'tab', 'data']
""" Available file formats for reading data. """


# --------------------------------------------------------------------------------
# Reading and Writing data files and result dictionaries
# --------------------------------------------------------------------------------


def write(filename_or_key, data, ext=None):
    """Write AnnData objects and dictionaries to file.

    If a key is passed, the filename is generated as
        filename = sett.writedir + key + sett.file_format_data
    This defaults to
        filename = 'write/' + key + '.h5'
    and can be changed by reseting sett.writedir and sett.file_format_data.

    Parameters
    ----------
    filename_or_key : str
        Filename of data file or key.
    data : dict, AnnData
        Annotated data object or dict storing arrays as values.
    """
    filename_or_key = str(filename_or_key)  # allow passing pathlib.Path objects
    from .classes.ann_data import AnnData
    if isinstance(data, AnnData):
        d = data.to_dict()
    else:
        d = data
    if is_filename(filename_or_key):
        filename = filename_or_key
        ext_ = is_filename(filename, return_ext=True)
        if ext is None:
            ext = ext_
        elif ext != ext_:
            raise ValueError('It suffices to provide the file type by '
                             'providing a proper extension to the filename.'
                             'One of "txt", "csv", "h5" or "npz".')
    else:
        key = filename_or_key
        filename = get_filename_from_key(key)
    ext = sett.file_format_data if ext is None else ext
    write_dict_to_file(filename, d, ext=ext)
    if ext in {'txt', 'csv'}:
        filename = get_filename_from_key(key, ext='h5')
        write_dict_to_file(filename, d, ext='h5')


def read(filename_or_key, sheet='', ext='', delim=None, first_column_names=None,
         as_strings=False, backup_url='', return_dict=False):
    """Read file or dictionary and return data dictionary.

    To speed up reading and save storage space, this creates an hdf5 file if
    it's not present yet.

    Parameters
    ----------
    filename_or_key : str, Path
        Filename of data file or key used in function write(key,dict).
    sheet : str, optional
        Name of sheet in Excel file.
    ext : str, optional (default: automatically inferred from filename)
        Extension that indicates the file type.
    delim : str, optional
        Delimiter that separates data within text file. If None, will split at
        arbitrary number of white spaces, which is different from enforcing
        splitting at single white space ' '.
    first_column_names : bool, optional
        Assume the first column stores samplenames. Is unnecessary if the sample
        names are not floats or integers: if they are strings, this will be
        detected automatically.
    as_strings : bool, optional
        Read names instead of numbers.
    backup_url : str, optional
        Retrieve the file from a URL if not present on disk.
    return_dict : bool, optional (default: False)
        Return dictionary instead of AnnData object.

    Returns
    -------
    data : sc.AnnData object or dict if return_dict == True
        If dict, it contains
        X : np.ndarray, optional
            Data array for further processing, columns correspond to genes,
            rows correspond to samples.
        row_names : np.ndarray, optional
            Array storing the names of rows (experimental labels of samples).
        col_names : np.ndarray, optional
            Array storing the names of columns (gene names).
        Maybe more items, if they are found in the file.
    """
    filename_or_key = str(filename_or_key)  # allow passing pathlib.Path objects
    from .classes.ann_data import AnnData
    if is_filename(filename_or_key):
        d = read_file(filename_or_key, sheet, ext, delim, first_column_names,
                      as_strings, backup_url)
        if return_dict:
            return d
        else:
            adata = AnnData(d)
            return adata

    # generate filename and read to dict
    key = filename_or_key
    filename = sett.writedir + key + '.' + sett.file_format_data
    if not os.path.exists(filename):
        raise ValueError('Reading with key ' + key + ' failed! ' +
                         'Provide valid key or valid filename directly: ' +
                         'inferred filename ' + filename + ' does not exist.\n' +
                         'If you intended to provide a filename, either ' +
                         'use a filename on one of the available extensions\n' +
                         str(avail_exts) +
                         '\nor provide the parameter "ext" to sc.read.')
    d = read_file_to_dict(filename, ext=sett.file_format_data)
    if return_dict:
        return d
    else:
        return AnnData(d)


# --------------------------------------------------------------------------------
# Reading and writing parameter files
# --------------------------------------------------------------------------------


def read_params(filename, asheader=False, verbosity=0):
    """
    Read parameter dictionary from text file.

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
    """
    Write parameters to file, so that it's readable by read_params.

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
    """
    Transform params list to dictionary.
    """
    params = {}
    for i in range(0, len(params_list)):
        key_val = params_list[i].split('=')
        if len(key_val) != 2:
            raise ValueError('Need to provide parameters as list of the form `par1=value1 par2=value2 ...`.')
        key, val = key_val
        params[key] = convert_string(val)
    return params


# --------------------------------------------------------------------------------
# Reading and Writing data files
# --------------------------------------------------------------------------------


def read_file(filename, sheet='', ext='', delim=None, first_column_names=None,
              as_strings=False, backup_url=''):
    """
    Read file and return data dictionary.

    To speed up reading and save storage space, this creates an hdf5 file if
    it's not present yet.

    Parameters
    ----------
    filename : str, Path
        Filename of data file.
    sheet : str, optional
        Name of sheet in Excel file.
    ext : str, optional (default: inferred from filename)
        Extension that indicates the file type.
    delim : str, optional
        Separator that separates data within text file. If None, will split at
        arbitrary number of white spaces, which is different from enforcing
        splitting at single white space ' '.
    first_column_names : bool, optional
        Assume the first column stores samplenames. Is unnecessary if the sample
        names are not floats or integers, but strings.
    as_strings : bool, optional
        Read names instead of numbers.
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
    if ext != '' and ext not in avail_exts:
        raise ValueError('Please provide one of the available extensions.\n'
                         + avail_exts)
    else:
        ext = is_filename(filename, return_ext=True)
    # check whether data file is present, otherwise download
    filename = check_datafile_present(filename, backup_url=backup_url)
    # read hdf5 files
    if ext == 'h5':
        if sheet == '':
            return read_file_to_dict(filename, ext=sett.file_format_data)
        else:
            sett.m(0, 'reading sheet', sheet, 'from file', filename)
            return _read_hdf5_single(filename, sheet)
    # read other file formats
    filename_stripped = filename.lstrip('./')
    if filename_stripped.startswith('data/'):
        filename_stripped = filename_stripped[5:]
    fast_ext = sett.file_format_data if sett.file_format_data in {'h5', 'npz'} else 'h5'
    filename_fast = (sett.writedir + 'data/'
                     + filename_stripped.replace('.' + ext, '.' + fast_ext))
    if not os.path.exists(filename_fast) or sett.recompute == 'read':
        sett.m(0, 'reading file', filename,
               '\n... writing an', sett.file_format_data,
               'version to speedup reading next time\n   ',
               filename_fast)
        if not os.path.exists(os.path.dirname(filename_fast)):
            os.makedirs(os.path.dirname(filename_fast))
        # do the actual reading
        if ext == 'xlsx' or ext == 'xls':
            if sheet == '':
                ddata = read_file_to_dict(filename, ext=ext)
            else:
                ddata = _read_excel(filename, sheet)
        elif ext == 'mtx':
            ddata = _read_mtx(filename)
        elif ext == 'csv':
            ddata = read_txt(filename, delim=',',
                             first_column_names=first_column_names,
                             as_strings=as_strings)
        elif ext in ['txt', 'tab', 'data']:
            if ext == 'data':
                sett.m(0, '... assuming ".data" means tab or white-space separated text file')
                sett.m(0, '--> change this by specifying ext to sc.read')
            ddata = read_txt(filename, delim, first_column_names,
                             as_strings=as_strings)
        elif ext == 'soft.gz':
            ddata = _read_softgz(filename)
        elif ext == 'txt.gz':
            sys.exit('TODO: implement similar to read_softgz')
        else:
            raise ValueError('Unkown extension', ext)
        # write as fast for faster reading when calling the next time
        write_dict_to_file(filename_fast, ddata, sett.file_format_data)
    else:
        ddata = read_file_to_dict(filename_fast, sett.file_format_data)
    return ddata

def _read_mtx(filename):
    """
    Read mtx file.
    """
    from scipy.io import mmread
    X = mmread(filename)
    from scipy.sparse.csr import csr_matrix
    X = csr_matrix(X)
    sett.m(0, '... did not find row_names or col_names')
    return {'X': X}

def read_txt(filename, delim=None, first_column_names=None, as_strings=False):
    """
    Return ddata dictionary.

    Parameters
    ----------
    filename : str, Path
        Filename to read from.
    delim : str, optional
        Separator that separates data within text file. If None, will split at
        arbitrary number of white spaces, which is different from enforcing
        splitting at single white space ' '.
    first_column_names : bool, optional
        Assume the first column stores samplenames.

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
    if as_strings:
        ddata = read_txt_as_strings(filename, delim)
    else:
        ddata = read_txt_as_floats(filename, delim, first_column_names)
    return ddata

def read_txt_as_floats(filename, delim=None, first_column_names=None):
    """
    Return data as list of lists of strings and the header as string.

    Parameters
    ----------
    filename : str, Path
        Filename of data file.
    delim : str, optional
        Separator that separates data within text file. If None, will split at
        arbitrary number of white spaces, which is different from enforcing
        splitting at single white space ' '.

    Returns
    -------
    data : list
         List of lists of strings.
    header : str
         String storing the comment lines (those that start with '#').
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
            line_list = line.split(delim)
            if not is_float(line_list[0]):
                col_names = line_list
                sett.m(0, '... assuming first line in file stores column names')
            else:
                if not is_float(line_list[0]) or first_column_names:
                    first_column_names = True
                    row_names.append(line_list[0])
                    data.append(line_list[1:])
                else:
                    data.append(line_list)
            break
    if not col_names:
        # try reading col_names from the last comment line
        if len(header) > 0:
            sett.m(0, '... assuming last comment line stores variable names')
            col_names = np.array(header.split('\n')[-2].strip('#').split())
        # just numbers as col_names
        else:
            sett.m(0, '... did not find column names in file')
            col_names = np.arange(len(data[0])).astype(str)
    col_names = np.array(col_names, dtype=str)
    # check if first column contains row names or not
    if first_column_names is None:
        first_column_names = False
    for line in f:
        line_list = line.split(delim)
        if not is_float(line_list[0]) or first_column_names:
            sett.m(0, '... assuming first column in file stores row names')
            first_column_names = True
            row_names.append(line_list[0])
            data.append(line_list[1:])
        else:
            data.append(line_list)
        break
    # parse the file
    for line in f:
        line_list = line.split(delim)
        if first_column_names:
            row_names.append(line_list[0])
            data.append(line_list[1:])
        else:
            data.append(line_list)
    sett.mt(0, 'read data into list of lists')
    # transfrom to array, this takes a long time and a lot of memory
    # but it's actually the same thing as np.genfromtext does
    # - we don't use the latter as it would involve another slicing step
    #   in the end, to separate row_names from float data, slicing takes
    #   a lot of memory and cpu time
    data = np.array(data, dtype=np.float64)
    sett.mt(0, 'constructed array from list of list')
    # transform row_names
    if not row_names:
        row_names = np.arange(len(data)).astype(str)
        sett.m(0, '... did not find row names in file')
    else:
        row_names = np.array(row_names)
        for iname, name in enumerate(row_names):
            row_names[iname] = name.strip('"')
    # adapt col_names if necessary
    if col_names.size > data.shape[1]:
        col_names = col_names[1:]
    for iname, name in enumerate(col_names):
        col_names[iname] = name.strip('"')
    ddata = {'X': data, 'row_names': row_names, 'col_names': col_names}
    return ddata

def read_txt_as_strings(filename, delim):
    """
    Interpret list of lists as strings
    """
    filename = str(filename)  # allow passing pathlib.Path objects
    # just a plain loop is enough
    header = ''
    data = []
    for line in open(filename):
        if line.startswith('#'):
            header += line
        else:
            line_list = line.split(delim)
            data.append(line_list)
    # now see whether we can simply transform it to an array
    if len(data[0]) == len(data[1]):
        X = np.array(data).astype(str)
        col_names = None
        row_names = None
        sett.m(0, '... the whole content of the file is in X')
    else:
        # strip quotation marks
        if data[0][0].startswith('"'):
            # using iterators over numpy arrays doesn't work efficiently here
            # speed no problem here, only done once
            for ir, r in enumerate(data):
                for ic, elem in enumerate(r):
                    data[ir][ic] = elem.strip('"')
        col_names = np.array(data[0]).astype(str)
        data = np.array(data[1:]).astype(str)
        row_names = data[:, 0]
        X = data[:, 1:]
        sett.m(0, '... first row is stored in "col_names"')
        sett.m(0, '... first column is stored in "row_names"')
        sett.m(0, '... data is stored in X')
    ddata = {'X': X, 'col_names': col_names, 'row_names': row_names}
    return ddata

def _read_hdf5_single(filename, key=''):
    """
    Read a single dataset from an hdf5 file.

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
        ddata = {'X' : X}
        # try to find row and column names
        for iname, name in enumerate(['row_names', 'col_names']):
            if name in keys:
                ddata[name] = f[name][()]
            elif key + '_' + name.replace('_', '') in keys:
                ddata[name] = f[key + '_' + name.replace('_', '')][()]
            else:
                ddata[name] = np.arange(X.shape[0 if name == 'row_names' else 1])
                if key == 'X':
                    sett.m(0, 'did not find', name, 'in', filename)
            ddata[name] = ddata[name].astype(str)
            if X.ndim == 1:
                break
    return ddata

def _read_excel(filename, sheet=''):
    """
    Read excel file and return data dictionary.

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
        df = read_excel(filename,sheet)
    except Exception as e:
        # in case this raises an error using Python 2.7
        print('if on Python 2.7 '
              'try installing xlrd using "sudo pip install xlrd"')
        raise e
    return ddata_from_df(df)


def _read_softgz(filename):
    """
    Read a SOFT format data file.

    The SOFT format is documented here
    http://www.ncbi.nlm.nih.gov/geo/info/soft2.html.

    Returns
    -------
    ddata : dict, containing
        X : np.ndarray
            A d x n array of gene expression values.
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
        I = [i for i,x in enumerate(sample_names) if x.startswith("GSM")]
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
            gene_names.append(#V[0] + ";" + # only use the second gene name
                              V[1])
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

# --------------------------------------------------------------------------------
# Reading and writing for dictionaries
# --------------------------------------------------------------------------------

def read_file_to_dict(filename, ext='h5'):
    """
    Read file and return dict with keys.

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
    sett.m(0, 'reading file', filename)
    d = {}
    if ext in {'h5', 'txt', 'csv'}:
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
    if 'X_csr_data' in d:
        d = load_sparse_csr(d)
    return d


def postprocess_reading(key, value):
    if value.dtype.kind == 'S':
        value = value.astype(str)
    return key, value


def preprocess_writing(key, value):
    if isinstance(key, dict):
        for k, v in value.items():
            return preprocess_writing(k, v)
    value = np.array(value)
    # some output about the data to write
    sett.m(1, key, type(value),
           value.dtype, value.dtype.kind, value.shape)
    # make sure string format is chosen correctly
    if value.dtype.kind == 'U':
        value = value.astype(np.string_)
    return key, value


def write_dict_to_file(filename, d, ext='h5'):
    """
    Write dictionary to file.

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
        sett.m(0, 'creating directory', directory + '/', 'for saving output files')
        os.makedirs(directory)
    if ext == 'h5' or ext == 'npz':
        d_write = {}
        from scipy.sparse import issparse
        for key, value in d.items():
            if key == 'X' and issparse(value):
                for k, v in save_sparse_csr(value).items():
                    d_write[k] = v
            else:
                key, value = preprocess_writing(key, value)
                d_write[key] = value
    if ext == 'h5':
        with h5py.File(filename, 'w') as f:
            for key, value in d_write.items():
                try:
                    f.create_dataset(key, data=value)
                except Exception as e:
                    sett.m(0, 'Error creating dataset for key =', key)
                    raise e
    elif ext == 'npz':
        np.savez(filename, **d_write)
    elif ext == 'csv' or ext == 'txt':
        # here this is actually a directory that corresponds to the
        # single hdf5 file
        dirname = filename.replace('.' + ext, '/')
        sett.m(0, '... exporting', ext, 'files to', dirname)
        if not os.path.exists(dirname):
            os.makedirs(dirname)
        for key, value in d.items():
            key, value = preprocess_writing(key, value)
            from pandas import DataFrame
            filename = dirname + '/' + key + '.' + ext
            if value.dtype.names is None:
                if value.dtype.char == 'S':
                    df = DataFrame(value.astype('U'))
                else:
                    df = DataFrame(value)
                df.to_csv(filename, sep=(' ' if ext == 'txt' else ','),
                          header=False, index=False)
            else:
                df = DataFrame.from_records(value)
                cols = list(df.select_dtypes(include=[object]).columns)
                # convert to unicode string
                df[cols] = df[cols].values.astype('U')
                if key == 'var':
                    df = df.T
                    df.to_csv(filename,
                              sep=(' ' if ext == 'txt' else ','),
                              header=False)
                else:
                    df.to_csv(filename,
                              sep=(' ' if ext == 'txt' else ','),
                              index=False)


# --------------------------------------------------------------------------------
# Type conversion
# --------------------------------------------------------------------------------


def save_sparse_csr(X):
    from scipy.sparse.csr import csr_matrix
    X = csr_matrix(X)
    return {'X_csr_data': X.data,
            'X_csr_indices': X.indices,
            'X_csr_indptr': X.indptr,
            'X_csr_shape': X.shape}


def load_sparse_csr(d):
    from scipy.sparse.csr import csr_matrix
    d['X'] = csr_matrix((d['X_csr_data'],
                         d['X_csr_indices'],
                         d['X_csr_indptr']),
                        shape=d['X_csr_shape'])
    del d['X_csr_data']
    del d['X_csr_indices']
    del d['X_csr_indptr']
    del d['X_csr_shape']
    return d


def is_float(string):
    """
    Check whether string is float.

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
    """
    Check whether string is integer.
    """
    try:
        int(string)
        return True
    except ValueError:
        return False


def convert_bool(string):
    """
    Check whether string is boolean.
    """
    if string == 'True':
        return True, True
    elif string == 'False':
        return True, False
    else:
        return False, False


def convert_string(string):
    """
    Convert string to int, float or bool.
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


# --------------------------------------------------------------------------------
# Helper functions for reading and writing
# --------------------------------------------------------------------------------


def get_filename_from_key(key, ext=None):
    ext = sett.file_format_data if ext is None else ext
    filename = sett.writedir + key + '.' + ext
    return filename


def ddata_from_df(df):
    """
    Write pandas.dataframe to ddata dictionary.
    """
    ddata = {
        'X': df.values[:, 1:].astype(float),
        'row_names': df.iloc[:, 0].values.astype(str),
        'col_names': np.array(df.columns[1:], dtype=str)}
    return ddata


def download_progress(count, blockSize, totalSize):
    percent = int(count*blockSize*100/totalSize)
    sys.stdout.write("\r" + "... %d%%" % percent)
    sys.stdout.flush()


def check_datafile_present(filename, backup_url=''):
    """
    Check whether the file is present, otherwise download.
    """
    if filename.startswith('sim/'):
        if not os.path.exists(filename):
            exkey = filename.split('/')[1]
            print('file ' + filename + ' does not exist')
            print('you can produce the datafile by')
            sys.exit('running subcommand "sim ' + exkey + '"')

    if not os.path.exists(filename):
        if os.path.exists('../' + filename):
            # we are in a subdirectory of the scanpy repo
            return '../' + filename
        else:
            # download the file
            if backup_url == '':
                sys.exit('file ' + filename + ' does not exist')
            sett.m(0, 'try downloading from url\n' + backup_url + '\n' +
                   '... this may take a while but only happens once')
            d = os.path.dirname(filename)
            if not os.path.exists(d):
                sett.m(0, 'creating directory', d+'/', 'for saving data')
                os.makedirs(d)
            from .compat.urllib_request import urlretrieve
            urlretrieve(backup_url, filename, reporthook=download_progress)
            sett.m(0, '')

    return filename


def is_filename(filename_or_key, return_ext=False):
    """ Check whether it is a filename. """
    for ext in avail_exts:
        l = len('.' + ext)
        # check whether it ends on the extension
        if '.'+ext in filename_or_key[-l:]:
            if return_ext:
                return ext
            else:
                return True
    if return_ext:
        raise ValueError('"' + filename_or_key + '"'
                         + ' does not contain a valid extension\n'
                         + 'Please provide one of the available extensions.\n'
                         + avail_exts)
    else:
        return False
