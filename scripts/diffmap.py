#!/usr/bin/env python
# coding: utf-8
"""
Wrapper for single Scanpy Tool: diffmap
=======================================

For package Scanpy (https://github.com/theislab/scanpy).
Written in Python 3 (compatible with 2).
"""   

from sys import path
# scanpy, first try loading it locally
path.insert(0, '.')
import scanpy as sc

sc.read_args_run_tool('diffmap')
