import sys
from setuptools import setup, find_packages
from distutils.extension import Extension
from pathlib import Path
import versioneer

# Pip calls setup.py egg_info to get static dependency information,
# then installs them, and finally it calls setup.py develop/bdist.
# We only need numpy when really installing/building the package,
# and we only rebuild the .pyx extension when cython is installed.
cmd_class = {}
ext_modules = []
include_dirs = []
if 'egg_info' not in sys.argv:
    import numpy

    include_dirs.append(numpy.get_include())

    try:
        from Cython.Distutils import build_ext
    except ImportError:
        suffix = 'c'
    else:
        cmd_class['build_ext'] = build_ext
        suffix = 'pyx'

    ext_modules.append(Extension(
        "scanpy.cython.utils_cy",
        ["scanpy/cython/utils_cy." + suffix],
        include_dirs=include_dirs,
    ))

package_name = 'scanpy'

req_path = Path('requires.txt')
if not req_path.is_file():
    req_path = Path(package_name + '.egg-info') / req_path
with req_path.open() as requirements:
    requires = [l.strip() for l in requirements]

with open('README.rst', encoding='utf-8') as readme_f:
    readme = readme_f.read()

setup(
    name=package_name,
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(cmd_class),
    description='Single-Cell Analysis in Python.',
    long_description=readme,
    url='http://github.com/theislab/scanpy',
    author='Alex Wolf, Philipp Angerer, Tobias Callies',
    author_email='alex.wolf@helmholtz-muenchen.de',
    license='BSD-3-Clause',
    entry_points={
        'console_scripts': [
            'scanpy = scanpy.__main__:main',
        ],
    },
    install_requires=requires,
    packages=find_packages(),  # + ['scanpy.sim_models'], might need to include sim_models
    include_dirs=include_dirs,
    # `package_data` does NOT work for source distributions!!!
    # you also need MANIFTEST.in
    # https://stackoverflow.com/questions/7522250/how-to-include-package-data-with-setuptools-distribute
    package_data={'': '*.txt'},
    include_package_data=True,
    ext_modules=ext_modules,
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
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Visualization',
    ],
)
