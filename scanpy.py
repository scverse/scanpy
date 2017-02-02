#!/usr/bin/env python
"""
Wrapper for scanpy.scanpy
=========================

For package Scanpy (https://github.com/theislab/scanpy).
Written in Python 3 (compatible with 2).   
"""
from sys import path

# scanpy, first try loading it locally
#path.insert(0, '.')
from scanpy.__main__ import main

main()
