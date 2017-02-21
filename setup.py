import sys
from setuptools import setup

more_requires = []
if sys.version_info[:2] < (3, 5):
    more_requires.append('configparser')

setup(
    name='scanpy',
    version='0.1',
    description='Single-Cell Analysis in Python.',
    url='http://github.com/theislab/scanpy',
    author='F. Alexander Wolf',
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
        'xlrd',  # for reading excel data
        'h5py',
        'scikit-learn',
    ] + more_requires,
    packages=['scanpy', 'scanpy.tools', 'scanpy.compat', 'scanpy.exs',
    'scanpy.preprocess', 'scanpy.classes'],
    zip_safe=False,
)
