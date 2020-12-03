from pathlib import Path

here = Path(__file__).parent

try:
    from setuptools_scm import get_version
    import pytoml

    proj = pytoml.loads((here.parent / 'pyproject.toml').read_text())
    metadata = proj['tool']['scanpy']

    __version__ = get_version(root='..', relative_to=__file__)
    __author__ = metadata['author']
    __email__ = metadata['author-email']
except (ImportError, LookupError, FileNotFoundError):
    from ._compat import pkg_metadata

    metadata = pkg_metadata(here.name)
    __version__ = metadata['Version']
    __author__ = metadata['Author']
    __email__ = metadata['Author-email']
