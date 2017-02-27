#!/usr/bin/env python
"""
Wrapper for scanpy.scanpy
"""
from sys import path
path.insert(0, '.')
path.insert(0, '..')
from scanpy.__main__ import main
main()
