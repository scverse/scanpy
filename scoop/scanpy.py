"""Add pli to scanpy"""

from . import pli

def patch():
    import scanpy as sc
    setattr(sc, 'pli', pli)

patch()