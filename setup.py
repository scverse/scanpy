import sys

if sys.version_info < (3, 6):
    sys.exit('scanpy requires Python >= 3.6')
from pathlib import Path

from setuptools import setup, find_packages

try:
    import pytoml
except ImportError:
    sys.exit('Please use `pip install .` or install pytoml first.')

proj = pytoml.loads(Path('pyproject.toml').read_text())
metadata = proj['tool']['scanpy']

setup(
    name='scanpy',
    use_scm_version=True,
    setup_requires=['setuptools_scm'],
    description='Single-Cell Analysis in Python.',
    long_description=Path('README.rst').read_text('utf-8'),
    url='http://github.com/theislab/scanpy',
    author=metadata['author'],
    author_email=metadata['author-email'],
    license='BSD',
    python_requires='>=3.6',
    install_requires=[
        l.strip() for l in Path('requirements.txt').read_text('utf-8').splitlines()
    ],
    extras_require=dict(
        louvain=['python-igraph', 'louvain>=0.6,!=0.6.2'],
        leiden=['python-igraph', 'leidenalg'],
        bbknn=['bbknn'],
        scvi=['scvi==0.6.7'],
        rapids=['cudf>=0.9', 'cuml>=0.9', 'cugraph>=0.9'],
        magic=['magic-impute>=2.0'],
        skmisc=['scikit-misc>=0.1.3'],
        harmony=['harmonypy'],
        scanorama=['scanorama'],
        scrublet=['scrublet'],
        dev=['setuptools_scm', 'pytoml', 'black>=20.8b1'],
        doc=[
            'sphinx>=3.2',
            'sphinx_rtd_theme>=0.3.1',
            'sphinx_autodoc_typehints',
            'scanpydoc>=0.5',
            'typing_extensions; python_version < "3.8"',  # for `Literal`
        ],
        test=[
            'pytest>=4.4',
            'dask[array]!=2.17.0',
            'fsspec',
            'zappy',
            'zarr',
            'profimp',
        ],
    ),
    packages=find_packages(),
    include_package_data=True,
    entry_points=dict(console_scripts=['scanpy=scanpy.cli:console_main']),
    zip_safe=False,
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Environment :: Console',
        'Framework :: Jupyter',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Natural Language :: English',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: Microsoft :: Windows',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Visualization',
    ],
)
