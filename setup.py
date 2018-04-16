from setuptools import setup, find_packages
from pathlib import Path
import versioneer

package_name = 'scanpy'

req_path = Path('requires.txt')
if not req_path.is_file():
    req_path = Path(package_name + '.egg-info') / req_path
with req_path.open() as requirements:
    requires = [l.strip() for l in requirements]

with open('README.rst', encoding='utf-8') as readme_f:
    readme = readme_f.read()

author = 'Alex Wolf, Philipp Angerer, Davide Cittaro, Gokcen Eraslan, Tobias Callies'

setup(
    name=package_name,
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description='Single-Cell Analysis in Python.',
    long_description=readme,
    url='http://github.com/theislab/scanpy',
    author=author,
    author_email='alex.wolf@helmholtz-muenchen.de',
    license='BSD',
    install_requires=requires,
    packages=find_packages(),
    # `package_data` does NOT work for source distributions!!!
    # you also need MANIFTEST.in
    # https://stackoverflow.com/questions/7522250/how-to-include-package-data-with-setuptools-distribute
    package_data={'': '*.txt'},
    include_package_data=True,
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
