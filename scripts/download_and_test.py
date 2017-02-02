#! /usr/bin/env python
"""
Download and Test Scanpy
========================

Run this script from any directory to download Scanpy from GitHub and run the
tests in test.sh.
"""
import os
import sys
import zipfile
import shutil

# download the master branch
if sys.version_info >= (3,0):
    from urllib.request import urlretrieve
else: # Python 2
    from urllib import urlretrieve

url = 'https://github.com/theislab/scanpy/archive/master.zip'
dirname = 'scanpy_test'
filename = dirname + '.zip'
print('downloading master branch from '+ url)
urlretrieve(url,filename)

print('unzipping ' + filename + ' to directory ' + dirname)
zip_ref = zipfile.ZipFile(filename, 'r')
dirname_tmp = 'test_scanpy_tmp'
zip_ref.extractall(dirname_tmp)
zip_ref.close()
shutil.move(dirname_tmp+'/scanpy-master',dirname)
os.rmdir(dirname_tmp)

os.chdir(dirname)
print('run tests in ./tests/test.sh')
for line in open('./tests/test.sh'):
    line = line.strip()
    if line.startswith('python'):
        width = 2 + len(line)
        print(width*'-')
        print('+ ' + line)
        print(width*'-')
        os.system(line)

# remove test directory
# os.chdir('..')
# os.remove(filename)
# shutil.rmtree(dirname)
