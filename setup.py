import sys
from setuptools import setup
from distutils.extension import Extension
import numpy

use_cython = False  # set this to False
if use_cython:
    from Cython.Distutils import build_ext

cmdclass = {}
ext_modules = []
if use_cython:
    ext_modules += [
        Extension("scanpy.cython.utils_cy",
                  ["scanpy/cython/utils_cy.pyx"]),
    ]
    cmdclass.update({ 'build_ext': build_ext })
else:
    ext_modules += [
        Extension("scanpy.cython.utils_cy",
                  ["scanpy/cython/utils_cy.c"]),
    ]

more_requires = []
# if sys.version_info[:2] < (3, 5):
more_requires.append('configparser')
# if sys.version_info[:2] < (3, 4):
more_requires.append('enum34')  # we specifically seem to need enum34
if sys.version_info[:2] < (3, 0):
    more_requires.append('xlrd')  # for reading excel data

setup(
    name='scanpy',
    version='0.1',
    description='Single-Cell Analysis in Python.',
    url='http://github.com/theislab/scanpy',
    author='F. Alexander Wolf, P. Angerer',
    author_email='alex.wolf@helmholtz-muenchen.de',
    license='GPL-3.0',
    entry_points={
        'console_scripts': [
            'scanpy = scanpy.__main__:main',
        ],
    },
    install_requires=[
        'matplotlib',
        'pandas',
        'scipy',
        'h5py',          # hdf5 file reading and writing
        'scikit-learn',  # standard machine-learning algorithms
        'statsmodels',   # standard statistical models
        'natsort',       # natural, human-readable sorting
        'joblib',        # simple parallel computing
        'profilehooks'   # profiling
    ] + more_requires,
    packages=[
        'scanpy',
        'scanpy.tools',
        'scanpy.compat',
        'scanpy.examples',
        'scanpy.preprocess',
        'scanpy.classes',
        'scanpy.sim_models',
        'scanpy.cython',
    ],
    include_dirs=[numpy.get_include()],
    cmdclass=cmdclass,
    ext_modules=ext_modules,
    zip_safe=False,
)
