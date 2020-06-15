import scipy.sparse as sp


def _check_rpy2():
    try:
        import rpy2
    except ImportError:
        raise ImportError("Please install rpy2 package.")


def _check_anndata2ri():
    try:
        import anndata2ri
    except ImportError:
        raise ImportError("Please install anndata2ri package.")


def _is_installed(package_name):
    _check_rpy2()

    from rpy2.robjects.packages import isinstalled

    if not isinstalled(package_name):
        raise ImportError(f"Please install {package_name} R package.")


def _set_seed(seed):
    _check_rpy2()

    from rpy2.robjects import r

    set_seed = r('set.seed')
    set_seed(seed)


def _py2r(x):
    _check_rpy2()
    _check_anndata2ri()

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


def _r2py(x):
    _check_rpy2()
    _check_anndata2ri()

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
