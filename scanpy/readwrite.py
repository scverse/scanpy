# Copyright 2016-2017 F. Alexander Wolf (http://falexwolf.de).
"""
Reading and Writing

TODO
----
- Preserve case when writing params to files.
- Consider using openpyxl or xlsxwriter instead of pandas for reading and
  writing Excel files.   
"""

import os
import sys
import h5py
import numpy as np

from . import settings as sett

avail_exts = ['csv','xlsx','txt','h5','soft.gz','txt.gz', 'mtx', 'tab', 'data']
""" Available file formats for reading data. """

#--------------------------------------------------------------------------------
# Reading and Writing data files and result dictionaries
#--------------------------------------------------------------------------------

def write(filename_or_key, dict_or_adata):
    """
    Writes dictionaries - as returned by tools - to file.
    
    If a key is specified, the filename is generated as
        filename = sett.writedir + key + sett.extd
    This defaults to
        filename = 'write/' + key + '.h5'
    and can be changed by reseting sett.writedir and sett.extd.

    Parameters
    ----------
    filename_or_key : str
        Filename of data file or key used in function write(key,dict).
    dict_or_adata : dict, AnnData
        Annotated data matrix or dictionary convertible to one.
    """
    from .ann_data import AnnData
    if isinstance(dict_or_adata, AnnData):
        dictionary = dict_or_adata.to_dict()
        dictionary['isadata'] = True
    else:
        dictionary = dict_or_adata
    if is_filename(filename_or_key):
        filename = filename_or_key
    else:
        key = filename_or_key
        filename = get_filename_from_key(key)
        if 'writekey' not in dictionary:
            dictionary['writekey'] = key
    write_dict_to_file(filename, dictionary, ext=sett.extd)

def read(filename_or_key, sheet='', ext='', delim=None, first_column_names=None,
         as_strings=False, backup_url=''):
    """
    Read file or dictionary and return data dictionary.

    To speed up reading and save storage space, this creates an hdf5 file if
    it's not present yet.

    Parameters
    ----------
    filename_or_key : str
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

    Returns
    -------
    data : AnnData, dict, if dict it usually contains
        X : np.ndarray, optional
            Data array for further processing, columns correspond to genes,
            rows correspond to samples.
        row_names : np.ndarray, optional
            Array storing the names of rows (experimental labels of samples).
        col_names : np.ndarray, optional
            Array storing the names of columns (gene names).
    """
    if is_filename(filename_or_key):
        return read_file(filename_or_key, sheet, ext, delim, first_column_names, 
                         as_strings, backup_url)

    # generate filename and read to dict
    key = filename_or_key
    filename = sett.writedir + key + '.' + sett.extd
    if not os.path.exists(filename):
        raise ValueError('Reading with key ' + key + ' failed! ' + 
                         'Provide valid key or valid filename directly: ' + 
                         'inferred filename ' +
                         filename + ' does not exist.\n' +
                         'If you intended to provide a filename, either ' + 
                         'use a filename on one of the available extensions\n' + 
                         str(avail_exts) + 
                         'or provide the parameter "ext" to sc.read.')
    d = read_file_to_dict(filename)
    if 'isadata' in d:
        from .ann_data import AnnData
        del d['isadata']
        return AnnData(d)
    else:
        return d

#--------------------------------------------------------------------------------
# Reading and writing parameter files
#--------------------------------------------------------------------------------

def read_params(filename, asheader=False):
    """ 
    Read parameter dictionary from text file.

    Assumes that parameters are specified in the format:  
        par1 = value1
        par2 = value2

    Comments that start with '#' are allowed.
    
    Parameters
    ----------
    filename : str
        Filename of data file.
    asheader : bool, optional
        Read the dictionary from the header (comment section) of a file.

    Returns
    -------
    params : dict
        Dictionary that stores parameters.
    """
    if not os.path.exists(filename):
        filename = '../' + filename
    if not asheader:
        sett.m(0,'reading params file',filename)
    from collections import OrderedDict as odict
    params = odict([])
    for line in open(filename):
        if '=' in line:
            if not asheader or line.startswith('#'):
                line = line[1:] if line.startswith('#') else line
                key, val = line.split('=')
                key = key.strip(); val = val.strip()
                params[key] = convert_string(val)
    return params

def write_params(filename, *args, **dicts):
    """
    Write parameters to file, so that it's readable py read_params.

    Uses INI file format.
    """
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
    if len(params_list)%2 != 0:
        raise ValueError('need to provide a list of key value pairs')
    params = {}
    for i in range(0,len(params_list),2):
        key, val = params_list[i:i+2]
        params[key] = convert_string(val)    
    return params

#--------------------------------------------------------------------------------
# Reading and Writing data files
#--------------------------------------------------------------------------------

def read_file(filename, sheet='', ext='', delim=None, first_column_names=None,
              as_strings=False, backup_url=''):
    """ 
    Read file and return data dictionary.

    To speed up reading and save storage space, this creates an hdf5 file if
    it's not present yet.

    Parameters
    ----------
    filename : str
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
    if ext != '': 
        if not ext in avail_exts:
            raise ValueError('Please provide one of the available extensions.'
                             + avail_exts)
    else:
        ext = is_filename(filename, return_ext=True)
    # check whether data file is present, otherwise download
    filename = check_datafile_present(filename, backup_url=backup_url)
    # actual reading
    if ext == 'h5':
        if sheet == '':
            d = read_file_to_dict(filename, ext='h5')
            if 'isadata' in d:
                from .ann_data import AnnData
                del d['isadata']
                return AnnData(d)
            else:
                return d
        sett.m(0, 'reading sheet', sheet, 'from file', filename)
        return _read_hdf5_single(filename, sheet)
    # if filename is not in the hdf5 format, do some more book keeping
    filename_hdf5 = filename.replace('.' + ext, '.h5')
    # just to make sure that we have a different filename
    if filename_hdf5 == filename:
        filename_hdf5 = filename + '.h5'
    if not os.path.exists(filename_hdf5) or sett.recompute == 'read':
        sett.m(0,'reading file', filename,
                 '\n--> write an hdf5 version to speedup reading next time')
        # do the actual reading
        if ext == 'xlsx' or ext == 'xls':
            if sheet=='':
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
                sett.m(0, 'assuming ".data" means tab or white-space separated text file')
                sett.m(0, '--> change this by specifying ext to sc.read')
            ddata = read_txt(filename, delim, first_column_names,
                               as_strings=as_strings)
        elif ext == 'soft.gz':
            ddata = _read_softgz(filename)
        elif ext == 'txt.gz':
            print('TODO: implement similar to read_softgz')
            sys.exit()
        else:
            raise ValueError('unkown extension', ext)

        # specify Scanpy type of dictionary
        if 'type' not in ddata:
            ddata['type'] = 'data'
        # write as hdf5 for faster reading when calling the next time
        write_dict_to_file(filename_hdf5, ddata)
    else:
        ddata = read_file_to_dict(filename_hdf5)

    return ddata

def _read_mtx(filename):
    """
    Read mtx file.
    """
    from scipy.io import mmread
    X = mmread(filename).todense()
    sett.m(0, '--> did not find any row_names or columnnames')
    sett.m(0, '--> TODO: don\'t use dense format')
    return {'X': X}

def read_txt(filename, delim=None, first_column_names=None, as_strings=False):
    """ 
    Return ddata dictionary.

    Parameters
    ----------
    filename : str
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
    filename : str
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
                sett.m(0, '--> assuming first line in file stores column names')
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
            sett.m(0,'--> assuming last comment line stores variable names')
            col_names = np.array(header.split('\n')[-2].strip('#').split())
        # just numbers as col_names
        else:
            sett.m(0,'--> did not find column names in file')
            col_names = np.arange(len(data[0])).astype(str)
    col_names = np.array(col_names, dtype=str)
    # check if first column contains row names or not
    if first_column_names is None:
        first_column_names = False
    for line in f:
        line_list = line.split(delim)
        if not is_float(line_list[0]) or first_column_names:
            sett.m(0, '--> assuming first column in file stores row names')
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
        sett.m(0,'--> did not find row names in file')
    else:
        row_names = np.array(row_names)
    # adapt col_names if necessary
    if col_names.size > data.shape[1]:
        col_names = col_names[1:]
    ddata = {'X': data, 'row_names': row_names, 'col_names': col_names}
    return ddata

def read_txt_as_strings(filename, delim):
    """
    Interpret list of lists as strings
    """
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
        sett.m(0,'--> the whole content of the file is in X')
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
        sett.m(0,'--> first row is stored in "col_names"')
        sett.m(0,'--> first column is stored in "row_names"')
        sett.m(0,'--> data is stored in X')
    ddata = {'X': X, 'col_names': col_names, 'row_names': row_names}
    return ddata

def _read_hdf5_single(filename, key=''):
    """ 
    Read a single dataset from an hdf5 file.

    See also function read_file_to_dict(), which reads all keys of the hdf5 file.

    Parameters
    ----------
    filename : str
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
        for iname, name in enumerate(['row_names','col_names']):
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
    filename : str
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
    ddata = {'X': X, 'row_names': row_names, 'col_names': col_names, 
             'row': {'groups': groups}}
    return ddata

#--------------------------------------------------------------------------------
# Reading and writing for dictionaries
#--------------------------------------------------------------------------------

def read_file_to_dict(filename, ext='h5'):
    """ 
    Read file and return dict with keys.
   
    The recommended format for this is hdf5.
 
    If reading from an Excel file, key names correspond to sheet names.

    Parameters
    ----------
    filename : str
        Filename of data file.
    ext : {'h5', 'xlsx'}, optional
        Choose file format. Excel is much slower.

    Returns
    -------
    d : dict
    """
    sett.m(0,'reading file',filename)
    d = {}
    if ext == 'h5':
        with h5py.File(filename, 'r') as f:
            for key in f.keys():
                # the '()' means 'read everything' (by contrast, ':' only works
                # if not reading a scalar type)
                value = f[key][()]
                if value.dtype.kind == 'S':
                    value = value.astype(str)
                # get back dictionaries
                if key.endswith('_dict'):
                    dd = {}
                    for row in value:
                        dd[row[0]] = row[1:]
                    d[key[:-5]] = dd
                else:
                    d[key] = value
    elif ext == 'xlsx':
        import pandas as pd
        xl = pd.ExcelFile(filename)
        for sheet in xl.sheet_names:
            d[sheet] = xl.parse(sheet).values
    return d

def write_dict_to_file(filename, d, ext='h5'):
    """ 
    Write content of dictionary to file.

    Parameters
    ----------
    filename : str
        Filename of data file.
    d : dict
        Dictionary storing keys with np.ndarray-like data or scalars.
    ext : string
        Determines file type, allowed are 'h5' (hdf5),
        'xlsx' (Excel) [or 'csv' (comma separated value file)].
    """
    directory = os.path.dirname(filename)
    if not os.path.exists(directory):
        os.makedirs(directory)
    if ext == 'h5':
        with h5py.File(filename, 'w') as f:
            for key, value in d.items():
                # if is a dict, build an array from it
                if isinstance(value, dict):
                    array = []
                    for k, v in value.items():
                        array.append(np.r_[np.array([k]), v])
                    value = np.array(array)
                    key = key + '_dict'
                if type(value) != np.ndarray:
                    value = np.array(value)
                # some output about the data to write
                sett.m(1, key, type(value), 
                       value.dtype, value.dtype.kind, value.shape)
                # make sure string format is chosen correctly
                if value.dtype.kind == 'U':
                    value = value.astype(np.string_)
                # try writing
                try:
                    f.create_dataset(key, data=value)
                except Exception as e:
                    sett.m(0,'error creating dataset for key =', key)
                    raise e
    elif ext == 'xlsx':
        import pandas as pd
        with pd.ExcelWriter(filename, engine='openpyxl') as writer:
            for key, value in d.items():
                pd.DataFrame(value).to_excel(writer,key)

#--------------------------------------------------------------------------------
# Type conversion
#--------------------------------------------------------------------------------

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
    else:
        return string

#--------------------------------------------------------------------------------
# Helper functions for reading and writing
#--------------------------------------------------------------------------------

def get_filename_from_key(key):
    filename = sett.writedir + key + '.' + sett.extd
    return filename

def ddata_from_df(df):
    """
    Write pandas.dataframe to ddata dictionary.
    """
    ddata = {
        'X': df.values[:,1:].astype(float),
        'row_names': df.iloc[:,0].values.astype(str),
        'col_names': np.array(df.columns[1:], dtype=str)
        }
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
            # we are in a subdirectory of Scanpy
            return '../' + filename
        else:
            # download the file
            sett.m(0, 'file ' + filename + ' is not present')
            if backup_url == '':
                sys.exit(0)
            sett.m(0, 'try downloading from url\n' + backup_url + '\n' + 
                  '... this may take a while but only happens once')
            d = os.path.dirname(filename)
            if not os.path.exists(d):
                os.makedirs(d)
            from .compat.urllib_request import urlretrieve
            urlretrieve(backup_url, filename, reporthook=download_progress)
            sett.m(0,'')

    return filename

def is_filename(filename_or_key,return_ext=False):
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
        raise ValueError(filename_or_key 
                         + ' does not contain a valid extension\n'
                         + 'choose a filename that contains one of\n'
                         + avail_exts)
    else:
        return False


