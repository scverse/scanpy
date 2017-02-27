#!/usr/bin/env python
"""
Wrapper for Scanpy tool
"""   
from sys import path
# scanpy, first try loading it locally
path.insert(0, '.')
path.insert(0, '..')
import scanpy as sc
sc._read_command_line_args_run_single_tool('diffmap')
