import functools
import scipy.sparse as sp


def rpy2_check(func):
    '''Decorator to check whether rpy2 is installed at runtime'''

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        try:
            import rpy2
        except ImportError:
            raise ImportError("Please install rpy2 package.")
        return func(*args, **kwargs)

    return wrapper


def anndata2ri_check(func):
    '''Decorator to check whether anndata2ri is installed at runtime'''

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        try:
            import anndata2ri
        except ImportError:
            raise ImportError("Please install anndata2ri package.")
        return func(*args, **kwargs)

    return wrapper


@rpy2_check
def r_is_installed(package_name):
    '''Checks whether a given R package is installed'''
    from rpy2.robjects.packages import isinstalled

    if not isinstalled(package_name):
        raise ImportError(f"Please install {package_name} R package.")


@rpy2_check
def r_set_seed(seed):
    '''Set the seed of R random number generator'''
    from rpy2.robjects import r

    set_seed = r('set.seed')
    set_seed(seed)


@rpy2_check
def r_set_logger_level(level):
    '''Set the logger level of rpy2'''
    import rpy2.rinterface_lib.callbacks

    rpy2.rinterface_lib.callbacks.logger.setLevel(level)


@rpy2_check
@anndata2ri_check
def py2r(x):
    '''Convert a Python object to an R object using rpy2'''
    import rpy2.robjects as ro
    from rpy2.robjects import numpy2ri, pandas2ri
    from rpy2.robjects.conversion import localconverter
    import anndata2ri

    if sp.issparse(x):
        # workaround for: https://github.com/theislab/anndata2ri/issues/47
        return anndata2ri.scipy2ri.py2rpy(x)

    with localconverter(
        ro.default_converter + numpy2ri.converter + pandas2ri.converter
    ):
        x = ro.conversion.py2rpy(x)

    return x


@rpy2_check
@anndata2ri_check
def r2py(x):
    '''Convert an rpy2 (R)  object to a Python object'''
    import rpy2.robjects as ro
    from rpy2.robjects import numpy2ri, pandas2ri
    from rpy2.robjects.conversion import localconverter
    import anndata2ri

    try:
        with localconverter(
            ro.default_converter
            + numpy2ri.converter
            + pandas2ri.converter
            + anndata2ri.scipy2ri.converter
        ):
            x = ro.conversion.rpy2py(x)

    except TypeError:
        # workaround for: https://github.com/theislab/anndata2ri/issues/47
        x = anndata2ri.scipy2ri.rpy2py(x)

    return x
