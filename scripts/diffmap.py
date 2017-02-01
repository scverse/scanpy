#!/usr/bin/env python
"""
Wrapper for Scanpy tool
"""   

from sys import path
# scanpy, first try loading it locally
path.insert(0, '.')
import scanpy as sc

sc.read_args_run_tool('diffmap')
