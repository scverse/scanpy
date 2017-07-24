from setuptools import setup, find_packages
from distutils.extension import Extension
import numpy
import versioneer

try:
    from Cython.Distutils import build_ext
except ImportError:
    use_cython = False
else:
    use_cython = True

cmdclass = {}
ext_modules = []
if use_cython:
    ext_modules += [
        Extension("scanpy.cython.utils_cy",
                  ["scanpy/cython/utils_cy.pyx"]),
    ]
    cmdclass.update({'build_ext': build_ext})
else:
    ext_modules += [
        Extension("scanpy.cython.utils_cy",
                  ["scanpy/cython/utils_cy.c"]),
]

with open('requirements.txt') as requirements:
    requires = [l.strip() for l in requirements]

with open('README.rst') as readme_f:
    readme = readme_f.read()

setup(
    name='scanpy',
    version=versioneer.get_version(),
    description='Single-Cell Analysis in Python.',
    long_description=readme,
    url='http://github.com/theislab/scanpy',
    author='F. Alexander Wolf, P. Angerer',
    author_email='alex.wolf@helmholtz-muenchen.de',
    license='BSD-3-Clause',
    entry_points={
        'console_scripts': [
            'scanpy = scanpy.__main__:main',
        ],
    },
    install_requires=requires,
    packages=find_packages(exclude=['scripts', 'scripts.*']),
    include_dirs=[numpy.get_include()],
    cmdclass=versioneer.get_cmdclass(cmdclass),
    ext_modules=ext_modules,
    zip_safe=False,
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Environment :: Console',
        'Framework :: Jupyter',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Natural Language :: English',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: Microsoft :: Windows',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Visualization',
    ],
)
