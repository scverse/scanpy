from setuptools import setup, find_packages
from pathlib import Path
import versioneer

package_name = 'scanpy'

req_path = Path('requirements.txt')
with req_path.open() as requirements:
    requires = [l.strip() for l in requirements]

print(requires)

with open('README.rst', encoding='utf-8') as readme_f:
    readme = readme_f.read()

author = 'Alex Wolf, Philipp Angerer, Davide Cittaro, Gokcen Eraslan, Fidel Ramirez, Tobias Callies'

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
    extras_require=dict(
        louvain=['python-igraph', 'louvain>=0.6'],
        doc=['sphinx', 'sphinx_rtd_theme', 'sphinx_autodoc_typehints'],
        test=['pytest'],
    ),
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
