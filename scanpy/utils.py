# Copyright 2016-2017 F. Alexander Wolf (http://falexwolf.de).
"""
Utility functions and classes
=============================

TODO
----
- Preserve case when writing params to files.
- Consider using openpyxl or xlsxwriter instead of pandas for reading and
  writing Excel files.   
"""   

# this is necessary to import scanpy from within package
from __future__ import absolute_import
# standard modules
import os
import argparse
import h5py
import sys
import warnings
import traceback
import gzip
from collections import OrderedDict as odict
# scientific modules
import numpy as np
import scipy as sp
import scipy.cluster
from .compat.matplotlib import pyplot as pl
# local modules
import scanpy as sc
from . import settings as sett

avail_exts = ['csv','xlsx','txt','h5','soft.gz','txt.gz']
""" Available file formats for writing data. """

#--------------------------------------------------------------------------------
# Deal with examples
#--------------------------------------------------------------------------------

def fill_in_datakeys(dexamples, dexdata):
    """
    Update the 'examples dictionary' _examples.dexamples.

    If a datakey (key in 'datafile dictionary') is not present in the 'examples
    dictionary' it is used to initialize an entry with that key.
    
    If not specified otherwise, any 'exkey' (key in 'examples dictionary') is
    used as 'datakey'.
    """
    # default initialization of 'datakey' key with entries from data dictionary
    for exkey in dexamples:
        if 'datakey' not in dexamples[exkey]:
            if exkey in dexdata:
                dexamples[exkey]['datakey'] = exkey
            else:
                dexamples[exkey]['datakey'] = 'unspecified in dexdata'
    return dexamples

#--------------------------------------------------------------------------------
# Deal with tool parameters
#--------------------------------------------------------------------------------

def init_params(params, default_params, check=True):
    """
    Update default paramaters with params.
    """
    _params = dict(default_params)
    if params: # allow for params to be None
        for key, val in params.items():
            if key in default_params:
                _params[key] = val
            elif check:
                raise ValueError('\'' + key 
                                 + '\' is not a valid parameter key, '
                                 + 'consider one of \n' 
                                 + str(list(default_params.keys())))
    return _params

#--------------------------------------------------------------------------------
# Command-line argument reading and processing
#--------------------------------------------------------------------------------

def add_args(p, dadd_args=None):
    """
    Add arguments to parser.

    Parameters
    -------
    dadd_args : dict
        Dictionary of additional arguments formatted as 
            {'arg': {'type': int, 'default': 0, ... }}
    """

    aa = p.add_argument_group('Tool parameters').add_argument
    aa('exkey',
       type=str, default='', metavar='exkey',
       help='Specify the "example key" (just a shorthand), which is used '
            'to look up a data dictionary and parameters. ' 
            'Use Scanpy subcommand "examples" to inspect possible values.')
    aa('plotkey',
       type=str, default='', metavar='plotkey', nargs='?',
       help='Specify plotting tool for visualization (default: tool dependent).')
    # example key default argument
    aa('-p', '--params',
       nargs='*', default=None, metavar='k v',
       help='Provide optional parameters as list, '
            'e.g., "sigma 5 knn True" for setting "sigma" and "knn". See possible '
            'keys in the function definition above (default: "").')
    # make sure there are is conflict with dadd_args
    if dadd_args is None or '--paramsfile' not in dadd_args:
        aa('--paramsfile',
           type=str, default='', metavar='pf',
           help='Alternatively, specify the path to a parameter file (default: "").')
    # arguments from dadd_args
    if dadd_args is not None:
        for key, val in dadd_args.items():
            if key != 'arg':
                aa(key, **val)

    aa = p.add_argument_group('Plotting').add_argument
    aa('-q', '--plotparams',
       nargs='*', default=None, metavar='k v',
       help='Provide specific plotting parameters as list, '
            'e.g., "layout 3d cmap viridis". ' 
            'See possible keys by calling "--plotparams help" (default: "").')

    aa = p.add_argument_group('Toolchain').add_argument
    aa('--pre',
       type=str, default='', metavar='tool',
       help='Tool whose output should be used as input, ' 
            'often a tool that detects subgroups (default: tool dependent).')

    # standard arguments
    p = sett.add_args(p)

    return p

def read_args_tool(toolkey, dexamples, tool_add_args=None):
    """
    Read args for single tool.
    """
    p = default_tool_argparser(sc.help(toolkey, string=True), dexamples)
    if tool_add_args is None:
        p = add_args(p)
    else:
        p = tool_add_args(p)
    args = vars(p.parse_args())
    args = sett.process_args(args)
    return args

def default_tool_argparser(description, dexamples):
    """
    Create default parser for single tools.
    """
    epilog = '\n'
    for k,v in sorted(dexamples.items()):
        epilog += '  ' + k + '\n'
    p = argparse.ArgumentParser(
        description=description,
        add_help=False,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=('available values for examples (exkey):'+epilog))
    return p

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
    params = odict([])
    for line in open(filename):
        if '=' in line:
            if not asheader or line.startswith('#'):
                line = line[1:] if line.startswith('#') else line
                key, val = line.split('=')
                key = key.strip(); val = val.strip()
                params[key] = convert_string(val)
    return params

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

def write_params(filename,*args,**dicts):
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

#--------------------------------------------------------------------------------
# Reading and Writing data files and dictionaries
#--------------------------------------------------------------------------------

def write(filename_or_key, dictionary):
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
    """
    if 'writekey' not in dictionary:
        dictionary['writekey'] = filename_or_key
    filename = sett.writedir + filename_or_key + '.' + sett.extd
    if is_filename(filename_or_key):
        filename = filename_or_key    
    write_dict_to_file(filename, dictionary)

def read(filename_or_key, sheet='', sep=None, first_column_names=False, 
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
    sep : str, optional
        Separator that separates data within text file. If None, will split at
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
    ddata : dict containing
        X : np.ndarray
            Data array for further processing, columns correspond to genes,
            rows correspond to samples.
        rownames : np.ndarray
            Array storing the names of rows (experimental labels of samples).
        colnames : np.ndarray
            Array storing the names of columns (gene names).
    """
    if is_filename(filename_or_key):
        return read_file(filename_or_key, sheet, sep, first_column_names, 
                         as_strings, backup_url)

    # generate filename and read to dict
    key = filename_or_key
    filename = sett.writedir + key + '.' + sett.extd
    if not os.path.exists(filename):
        raise ValueError('Reading with key ' + key + ' failed! ' + 
                         'Provide valid key or filename directly: ' + 
                         'inferred filename ' +
                         filename + ' does not exist.')
    return read_file_to_dict(filename)

#--------------------------------------------------------------------------------
# Reading and Writing data files
#--------------------------------------------------------------------------------

def read_file(filename, sheet='', sep=None, first_column_names=False,
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
    sep : str, optional
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
        rownames : np.ndarray
            Array storing the names of rows (experimental labels of samples).
        colnames : np.ndarray
            Array storing the names of columns (gene names).

    If sheet is unspecified, and an h5 or xlsx file is read, the dict
    contains all sheets instead.
    """
    ext = is_filename(filename, return_ext=True)
    filename = check_datafile_present(filename, backup_url=backup_url)
    return _read_file(filename, ext, sheet, sep, first_column_names, 
                      as_strings)

def _read_file(filename, ext, sheet, sep, first_column_names,
               as_strings):
    """
    Same as read_file(), one additional parameter ext.

    Parameters
    ----------
    ext : str
        Extension that denotes the type of the file.
    """
    if ext == 'h5':
        if sheet == '':
            return read_file_to_dict(filename, ext='h5')
        sett.m(0, 'reading sheet', sheet, 'from file', filename)
        return _read_hdf5_single(filename, sheet)
    # if filename is not in the hdf5 format, do some more book keeping
    filename_hdf5 = filename.replace('.' + ext, '.h5')
    # just to make sure that we have a different filename
    if filename_hdf5 == filename:
        filename_hdf5 = filename + '.h5'
    if not os.path.exists(filename_hdf5):
        sett.m(0,'reading file', filename,
                 '\n--> write an hdf5 version to speedup reading next time')
        # do the actual reading
        if ext == 'xlsx' or ext == 'xls':
            if sheet=='':
                ddata = read_file_to_dict(filename, ext=ext)
            else:
                ddata = _read_excel(filename, sheet)
        elif ext == 'csv':
            ddata = _read_text(filename, sep=',', 
                               first_column_names=first_column_names,
                               as_strings=as_strings)
        elif ext == 'txt':
            ddata = _read_text(filename, sep, first_column_names,
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
        write_dict_to_file(filename_hdf5,ddata)
    else:
        ddata = read_file_to_dict(filename_hdf5)

    return ddata

def _read_text(filename, sep=None, first_column_names=False, as_strings=False):
    """ 
    Return ddata dictionary.

    Parameters
    ----------
    filename : str
        Filename to read from.
    sep : str, optional
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
        rownames : np.ndarray
            Array storing the names of rows (experimental labels of samples).
        colnames : np.ndarray
            Array storing the names of columns (gene names).
    """
    data, header = _read_text_raw(filename,sep)
    
    if as_strings:
        return _interpret_as_strings(data)
    else:
        return _interpret_as_floats(data, header, first_column_names)

def _read_text_raw(filename, sep=None):
    """
    Return data as list of lists of strings and the header as string.

    Parameters
    ----------
    filename : str
        Filename of data file.
    sep : str, optional
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
    for line in open(filename):
        if line.startswith('#'):
            header += line
        else:
            line_list = line.split(sep)
            data.append(line_list)
    return data, header

def _interpret_as_strings(data):
    """
    Interpret list of lists as floats.
    """
    if len(data[0]) == len(data[1]):
        X = np.array(data).astype(str)
    else:
        # strip quotation marks
        if data[0][0].startswith('"'):
            # using iterators over numpy arrays doesn't work efficiently here
            # speed no problem here, only done once
            for ir, r in enumerate(data):
                for ic, elem in enumerate(r):
                    data[ir][ic] = elem.strip('"')
        colnames = np.array(data[0]).astype(str)
        data = np.array(data[1:]).astype(str)
        rownames = data[:, 0]
        X = data[:, 1:]
    ddata = {'X': X, 'colnames': colnames, 'rownames': rownames}
    return ddata

def _interpret_as_floats(data, header, first_column_names):
    """
    Interpret as float array with optional colnames and rownames.
    """
    # if the first element of the data list cannot be interpreted as float, the
    # first row of the data is assumed to store variable names
    if not is_float(data[0][0]):
        sett.m(0,'--> assuming first line in file stores variable names')
        # if the first row is one element shorter
        colnames = np.array(data[0]).astype(str)
        data = np.array(data[1:])
    # try reading colnames from the last comment line
    elif len(header) > 0 and type(header) == str:
        sett.m(0,'--> assuming last comment line stores variable names')
        potentialnames = header.split('\n')[-2].strip('#').split()
        colnames = np.array(potentialnames)
        # skip the first column
        colnames = colnames[1:]
        data = np.array(data)
    # just numbers as colnames
    else:
        sett.m(0,'--> did not find variable names in file')
        data = np.array(data)
        colnames = np.arange(data.shape[1]).astype(str)

    # if the first element of the second row of the data cannot be interpreted
    # as float, it is assumed to store sample names
    if not is_float(data[1][0]) or first_column_names: 
        sett.m(0,'--> assuming first column stores sample names')
        rownames = data[:,0].astype(str)
        # skip the first column
        try:
            X = data[:, 1:].astype(float)
        except ValueError as e:
            msg = 'formating is strange, last line of data matrix is \n'
            msg += str(data[-1, 1:])
            raise ValueError(msg)
        if colnames.size > X.shape[1]:
            colnames = colnames[1:]
    # just numbers as rownames
    else:
        sett.m(0,'--> did not find sample names in file')
        X = data.astype(float)
        rownames = np.arange(X.shape[0]).astype(str)

    ddata = {
        'X' : X, 'rownames' : rownames, 'colnames' : colnames
    }

    return ddata

def _read_hdf5_single(filename,key=''):
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
        rownames : np.ndarray
            Array storing the names of rows (experimental labels of samples).
        colnames : np.ndarray
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
        ddata = { 'X' : X }
        # set row and column names
        for iname, name in enumerate(['rownames','colnames']):
            if name in keys:
                ddata[name] = f[name][()]
            elif key + '_' + name in keys:
                ddata[name] = f[key + '_' + name][()]
            else:
                ddata[name] = np.arange(X.shape[0 if name == 'rownames' else 1])
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
        rownames : np.ndarray
            Array storing the names of rows (experimental labels of samples).
        colnames : np.ndarray
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

    The following values can be exported:
        GID : A list of gene identifiers of length d.
        SID : A list of sample identifiers of length n.
        STP : A list of sample desriptions of length d.
        X   : A dxn array of gene expression values.

    We translate this to the conventions of scanpy.

    Note
    ----
    The function is based on a script by Kerby Shedden.
    http://dept.stat.lsa.umich.edu/~kshedden/Python-Workshop/gene_expression_comparison.html
    """

    with gzip.open(filename) as fid:

        #print(help(fid))

        # The header part of the file contains information about the
        # samples. Read that information first.
        SIF = {}
        for line in fid:
            line = line.decode("utf-8")
            if line.startswith("!dataset_table_begin"):
                break
            elif line.startswith("!subset_description"):
                subset_description = line.split("=")[1].strip()
            elif line.startswith("!subset_sample_id"):
                subset_ids = line.split("=")[1].split(",")
                subset_ids = [x.strip() for x in subset_ids]
                for k in subset_ids:
                    SIF[k] = subset_description

        # Next line is the column headers (sample id's)
        SID = fid.readline().decode("utf-8").split("\t")

        # The column indices that contain gene expression data
        I = [i for i,x in enumerate(SID) if x.startswith("GSM")]

        # Restrict the column headers to those that we keep
        SID = [SID[i] for i in I]

        # Get a list of sample labels
        STP = [SIF[k] for k in SID]

        # Read the gene expression data as a list of lists, also get the gene
        # identifiers
        GID,X = [],[]
        for line in fid:
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
            GID.append(#V[0] + ";" + 
                       V[1])

    # Convert the Python list of lists to a Numpy array and transpose to match
    # the Scanpy convention of storing observations in rows and variables in
    # colums.
    X = np.array(X).T
    # rownames are the sample identifiers
    rownames = SID
    # labels identifying sets
    setlabels = STP
    # column names are the gene identifiers
    colnames = GID

    ddata = {
        'X' : X, 'rownames' : rownames, 'colnames' : colnames,
        'groupnames_n' : setlabels
    }

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
        Returns OrderedDict.
    """
    sett.m(0,'reading file',filename)
    d = odict([])
    if ext == 'h5':
        with h5py.File(filename, 'r') as f:
            for key in f.keys():
                # the '()' means 'read everything' (by contrast, ':' only works
                # if not reading a scalar type)
                value = f[key][()]
                if value.dtype.kind == 'S':
                    d[key] = value.astype(str)
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
                if type(value) != np.ndarray:
                    value = np.array(value)
                # some output about the data to write
                sett.m(1,key,type(value),value.dtype,value.dtype.kind,value.shape)
                # make sure string format is chosen correctly
                if value.dtype.kind == 'U':
                    value = value.astype(np.string_)
                # try writing
                try:
                    f.create_dataset(key,data=value)
                except Exception as e:
                    sett.m(0,'error creating dataset for key =', key)
                    raise e
    elif ext == 'xlsx':
        import pandas as pd
        with pd.ExcelWriter(filename,engine='openpyxl') as writer:
            for key, value in d.items():
                pd.DataFrame(value).to_excel(writer,key)

#--------------------------------------------------------------------------------
# Helper functions for reading and writing
#--------------------------------------------------------------------------------

def ddata_from_df(df):
    """
    Write pandas.dataframe to ddata dictionary.
    """
    ddata = { 
        'X' : df.values[:,1:].astype(float), 
        'rownames' : df.iloc[:,0].values.astype(str), 
        # TODO: check whether we always have to start with the 1st
        #       column, same as in csv file
        'colnames' : np.array(df.columns[1:],dtype=str) 
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
            exit('running subcommand "sim ' + exkey + '"')

    if not os.path.exists(filename):
        if os.path.exists('../' + filename):
            # we are in a subdirectory of Scanpy
            return '../' + filename
        else:
            # download the file
            sett.m(0,'file ' + filename + ' is not present\n' + 
                  'try downloading from url\n' + backup_url + '\n' + 
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

#--------------------------------------------------------------------------------
# Others
#--------------------------------------------------------------------------------

def pretty_dict_string(d, indent=0):
    """
    Pretty output of nested dictionaries.
    """
    s = ''
    for key, value in sorted(d.items()):
        s += '    ' * indent + str(key)
        if isinstance(value, dict):
             s += '\n' + pretty_dict_string(value, indent+1)
        else:
             s += ' = ' + str(value) + '\n'
    return s

def markdown_dict_string(d):
    """
    Markdown output that can be pasted in the examples/README.md.
    """
    # sort experimental data from simulated data
    sorted_keys = []
    sim_keys = []
    for key, value in sorted(d.items()):
        if 'type' in value:
            if 'sim' in value['type']:
                sim_keys.append(key)
        else:
            sorted_keys.append(key)
    len_exp = len(sorted_keys) - 1
    sorted_keys += sim_keys
    # format output
    s = 'Examples using experimental data.\n'
    for ikey, key in enumerate(sorted_keys):
        value = d[key]
        s += '* [' + key + '](#' + key + ')'
        if 'ref' in value:
            if 'doi' in value:
                link = 'http://dx.doi.org/' + value['doi']
            elif 'url' in value:
                link = value['url']
            s += (' - [' +  value['ref'].replace('et al.','*et al.*') 
                         + '](' + link +  ')')
        if 'title' in value:
            s += '   \n*' + value['title'] + '*'
        s += '\n'
        if ikey == len_exp:
            s += '\nExamples using simulated data.\n'
    return s

def merge_dicts(*dicts):
    """
    Given any number of dicts, shallow copy and merge into a new dict,
    precedence goes to key value pairs in latter dicts.

    Note
    ----
    http://stackoverflow.com/questions/38987/how-to-merge-two-python-dictionaries-in-a-single-expression
    """
    result = {}
    for d in dicts:
        result.update(d)
    return result

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

def masks(list_of_index_lists,n):
    """
    Make an array in which rows store 1d mask arrays from list of index lists.

    Parameters
    ----------
    n : int
        Maximal index / number of samples.
    """
    # make a list of mask arrays, it's easier to store 
    # as there is a hdf5 equivalent
    for il,l in enumerate(list_of_index_lists):
        mask = np.zeros(n,dtype=bool)
        mask[l] = True
        list_of_index_lists[il] = mask
    # convert to arrays
    masks = np.array(list_of_index_lists)
    return masks

def warn_with_traceback(message, category, filename, lineno, file=None, line=None):
    """
    Get full tracebacks when warning is raised by setting

    warnings.showwarning = warn_with_traceback
    
    See also
    --------
    http://stackoverflow.com/questions/22373927/get-traceback-of-warnings
    """
    traceback.print_stack()
    log = file if hasattr(file,'write') else sys.stderr
    sett.write(warnings.formatwarning(message, category, filename, lineno, line))

def transpose_ddata(ddata):
    """
    Transpose a data dictionary.

    Parameters
    ----------
    ddata : dict containing (at least)
        X : np.ndarray
            Data array for further processing, columns correspond to genes,
            rows correspond to samples.
        rownames : np.ndarray
            Array storing the names of rows.
        colnames : np.ndarray
            Array storing the names of columns.
    Returns
    -------
    ddata : dict 
        With X transposed and rownames and colnames interchanged.
    """
    ddata['X'] = ddata['X'].T
    colnames = ddata['colnames']
    ddata['colnames'] = ddata['rownames']
    ddata['rownames'] = colnames
    return ddata

def subsample(X,subsample=1,seed=0):
    """ 
    Subsample a fraction of 1/subsample samples from the rows of X.

    Parameters
    ----------
    X : np.ndarray
        Data array.
    subsample : int
        1/subsample is the fraction of data sampled, n = X.shape[0]/subsample.
    seed : int
        Seed for sampling.

    Returns
    -------
    Xsampled : np.ndarray
        Subsampled X.
    rows : np.ndarray 
        Indices of rows that are stored in Xsampled.
    """
    if subsample == 1 and seed == 0:
        return X, np.arange(X.shape[0],dtype=int)
    if seed == 0:
        # this sequence is defined simply by skipping rows
        # is faster than sampling
        rows = np.arange(0,X.shape[0],subsample,dtype=int)
        n = rows.size
        Xsampled = np.array(X[rows])
    if seed > 0:
        n = int(X.shape[0]/subsample)
        np.random.seed(seed)
        Xsampled, rows = subsample_n(X,n=n)
    sett.m(0,'subsampled to',n,'of',X.shape[0],'data points')
    return Xsampled, rows

def subsample_n(X,n=0,seed=0):
    """ 
    Subsample n samples from rows of array.

    Parameters
    ----------
    X : np.ndarray
        Data array.
    seed : int
        Seed for sampling.

    Returns
    -------
    Xsampled : np.ndarray
        Subsampled X.
    rows : np.ndarray 
        Indices of rows that are stored in Xsampled.
    """
    if n < 0:
        raise ValueError('n must be greater 0')
    np.random.seed(seed)
    n = X.shape[0] if (n == 0 or n > X.shape[0]) else n
    rows = np.random.choice(X.shape[0],size=n,replace=False)
    Xsampled = np.array(X[rows])
    return Xsampled, rows

def comp_distance(X,metric='euclidean'):
    """ 
    Compute distance matrix for data array X
    
    Parameters
    ----------
    X : np.ndarray
        Data array (rows store samples, columns store variables).
    metric : string
        For example 'euclidean', 'sqeuclidean', see sp.spatial.distance.pdist.

    Returns
    -------
    D : np.ndarray
        Distance matrix.
    """
    D = sp.spatial.distance.pdist(X,metric=metric)
    D = sp.spatial.distance.squareform(D)
    sett.mt(0,'computed distance matrix with metric =', metric)
    sett.m(5,D)
    if False:
        pl.matshow(D)
        pl.colorbar()
    return D

def hierarch_cluster(M):
    """ 
    Cluster matrix using hierarchical clustering.

    Parameters
    ----------
    M : np.ndarray
        Matrix, for example, distance matrix.

    Returns
    -------
    Mclus : np.ndarray 
        Clustered matrix.
    indices : np.ndarray
        Indices used to cluster the matrix.
    """
    link = sp.cluster.hierarchy.linkage(M)
    indices = sp.cluster.hierarchy.leaves_list(link)
    Mclus = np.array(M[:,indices])
    Mclus = Mclus[indices,:]
    sett.mt(0,'clustered matrix')
    if False:
        pl.matshow(Mclus)
        pl.colorbar()
    return Mclus, indices

def check_datafile_deprecated(filename, ext=None):
    """
    Check whether the file is present and is not just a placeholder.

    If the file is not present at all, look for a file with the same name but
    ending on '_url.txt', and try to download the datafile from there.

    If the file size is below 500 Bytes, assume the file is just a placeholder
    for a link to github. Download from github in this case.
    """
    if filename.startswith('sim/'):
        if not os.path.exists(filename):
            exkey = filename.split('/')[1]
            print('file ' + filename + ' does not exist')
            print('you can produce the datafile by')
            exit('running subcommand "sim ' + exkey + '"')
    if ext is None:
        _, ext = os.path.splitext()
    if ext.startswith('.'):
        ext = ext[1:]
    if not os.path.exists(filename) and not os.path.exists('../'+filename):
        basename = filename.replace('.'+ext,'')
        urlfilename = basename + '_url.txt'
        if not os.path.exists(urlfilename):
            urlfilename = '../' + urlfilename
        if not os.path.exists(urlfilename):
            raise ValueError('Neither ' + filename + 
                             ' nor ../' + filename +
                             ' nor files with a url to download from exist \n' + 
                             '--> move your data file to one of these places \n' +
                             '    or cd/browse into scanpy root directory.')
        with open(urlfilename) as f:
            url = f.readline().strip()
        sett.m(0,'data file is not present \n' + 
              'try downloading data file from url \n' + url + '\n' + 
              '... this may take a while but only happens once')
        # download the file
        urlretrieve(url,filename,reporthook=download_progress)
        sett.m(0,'')

    rel_filename = filename
    if not os.path.exists(filename):
        rel_filename = '../' + filename

    # if file is smaller than 500 Bytes = 0.5 KB
    threshold = 500
    if os.path.getsize(rel_filename) < threshold:
        # download the file
        # note that this has 'raw' in the address
        github_baseurl = r'https://github.com/theislab/scanpy/raw/master/'
        fileurl = github_baseurl + filename
        sett.m(0,'size of file',rel_filename,'is below',threshold/1000.,' kilobytes') 
        sett.m(0,'--> presumably is a placeholder for a git-lfs file')
        sett.m(0,'... if you installed git-lfs, you can use \'git lfs checkout\'')
        sett.m(0,'... to update all files with their actual content') 
        sett.m(0,'--> downloading data file from github using url')
        sett.m(0,fileurl)
        sett.m(0,'... this may take a while but only happens once')
        # make a backup of the small file
        from shutil import move
        basename = rel_filename.replace('.'+ext,'')
        move(rel_filename,basename+'_placeholder.txt')
        # download the file
        try:
            urlretrieve(fileurl,rel_filename,reporthook=download_progress)
            sett.m(0,'')
        except Exception as e:
            sett.m(0,e)
            sett.m(0,'when calling urlretrieve() in module scanpy.utils.utils')
            sett.m(0,'--> is the github repo/url private?')
            sett.m(0,'--> if you have access, manually download the file')
            sett.m(0,fileurl)
            sett.m(0,'replace the small file',rel_filename,'with the downloaded file')
            quit()

    return rel_filename
