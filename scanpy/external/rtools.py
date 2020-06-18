import functools
import scipy.sparse as sp


def rpy2_import(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        try:
            import rpy2
        except ImportError:
            raise ImportError("Please install rpy2 package.")
        return func(*args, **kwargs)

    return wrapper


def anndata2ri_import(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        try:
            import anndata2ri
        except ImportError:
            raise ImportError("Please install anndata2ri package.")
        return func(*args, **kwargs)

    return wrapper


@rpy2_import
def _is_installed(package_name):
    from rpy2.robjects.packages import isinstalled

    if not isinstalled(package_name):
        raise ImportError(f"Please install {package_name} R package.")


@rpy2_import
def _set_seed(seed):
    from rpy2.robjects import r

    set_seed = r('set.seed')
    set_seed(seed)


@rpy2_import
def _set_logger_level(level):
    import rpy2.rinterface_lib.callbacks

    rpy2.rinterface_lib.callbacks.logger.setLevel(level)


@rpy2_import
@anndata2ri_import
def _py2r(x):
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


@rpy2_import
@anndata2ri_import
def _r2py(x):
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
