import sys
from setuptools import setup
from distutils.extension import Extension
import numpy
import versioneer
from Cython.Distutils import build_ext

with open('requirements.txt') as requirements:
    requires = [l.strip() for l in requirements]

more_requires = []
if sys.version_info[0] == 2:
    more_requires = [
        'configparser',  # named ConfigParser in py2
        'xlrd',          # pandas on py2 reads XSL files with that
        'enum34',        # enum module introduced in python 3.4
    ]

setup(
    name='scanpy',
    version=versioneer.get_version(),
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
    install_requires=requires + more_requires,
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
    # cmdclass={
    #     'build_ext': build_ext,
    # },
    cmdclass=versioneer.get_cmdclass({'build_ext': build_ext}),
    ext_modules=[
        Extension("scanpy.cython.utils_cy",
                  ["scanpy/cython/utils_cy.pyx"]),
    ],
    zip_safe=False,
)
