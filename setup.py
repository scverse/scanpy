import sys
<<<<<<< HEAD
if sys.version_info < (3,):
    sys.exit('scanpy requires Python >= 3.5')
from pathlib import Path

from setuptools import setup, find_packages
import versioneer
=======

if sys.version_info < (3, 6):
    sys.exit('scanpy requires Python >= 3.6')
from pathlib import Path

from setuptools import setup, find_packages
>>>>>>> upstream/master


try:
    from scanpy import __author__, __email__
except ImportError:  # Deps not yet installed
    __author__ = __email__ = ''

setup(
    name='scanpy',
<<<<<<< HEAD
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
=======
    use_scm_version=True,
    setup_requires=['setuptools_scm'],
>>>>>>> upstream/master
    description='Single-Cell Analysis in Python.',
    long_description=Path('README.rst').read_text('utf-8'),
    url='http://github.com/theislab/scanpy',
    author=__author__,
    author_email=__email__,
    license='BSD',
<<<<<<< HEAD
    python_requires='>=3.5',
    install_requires=[
        l.strip() for l in
        Path('requirements.txt').read_text('utf-8').splitlines()
=======
    python_requires='>=3.6',
    install_requires=[
        l.strip()
        for l in Path('requirements.txt').read_text('utf-8').splitlines()
>>>>>>> upstream/master
    ],
    extras_require=dict(
        louvain=['python-igraph', 'louvain>=0.6'],
        leiden=['python-igraph', 'leidenalg'],
        bbknn=['bbknn'],
<<<<<<< HEAD
        doc=['sphinx', 'sphinx_rtd_theme', 'sphinx_autodoc_typehints', 'scanpydoc'],
        test=['pytest>=3.9'],
    ),
    packages=find_packages(),
    # `package_data` does NOT work for source distributions!!!
    # you also need MANIFTEST.in
    # https://stackoverflow.com/questions/7522250/how-to-include-package-data-with-setuptools-distribute
    package_data={'': '*.txt'},
    include_package_data=True,
=======
        rapids=['cudf', 'cuml', 'cugraph'],
        doc=[
            'sphinx',
            'sphinx_rtd_theme',
            'sphinx_autodoc_typehints',
            'scanpydoc>=0.4.3',
            'typing_extensions; python_version < "3.8"',  # for `Literal`
        ],
        test=[
            'pytest>=4.4',
            'dask[array]',
            'fsspec',
            'zappy',
            'zarr',
            'black',
            'profimp',
        ],
    ),
    packages=find_packages(),
    entry_points=dict(console_scripts=['scanpy=scanpy.cli:console_main']),
>>>>>>> upstream/master
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
