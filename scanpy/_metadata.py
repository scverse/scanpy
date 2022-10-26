import traceback
from pathlib import Path

here = Path(__file__).parent


def refresh_entry_points():
    """\
    Under some circumstances, (e.g. when installing a PEP 517 package via pip),
    pkg_resources.working_set.entries is stale. This tries to fix that.
    See https://github.com/pypa/setuptools_scm/issues/513
    """
    try:
        import sys
        import pkg_resources

        ws: pkg_resources.WorkingSet = pkg_resources.working_set
        for entry in sys.path:
            ws.add_entry(entry)
    except Exception:
        pass


try:
    from setuptools_scm import get_version

    refresh_entry_points()
    __version__ = get_version(root='..', relative_to=__file__)
except (ImportError, LookupError, FileNotFoundError):
    from ._compat import pkg_metadata

    metadata = pkg_metadata(here.name)
    __version__ = metadata['Version']


def within_flit():
    """\
    Checks if we are being imported by flit.
    This is necessary so flit can import __version__ without all depedencies installed.
    There are a few options to make this hack unnecessary, see:
    https://github.com/takluyver/flit/issues/253#issuecomment-737870438
    """
    for frame in traceback.extract_stack():
        if frame.name == 'get_docstring_and_version_via_import':
            return True
    return False
