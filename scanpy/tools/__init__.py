"""
Tools
"""

from . import diffmap
from . import difftest
from . import dpt
from . import tsne
from . import sim
from . import drawg

try:
    # development tools
    from . import ctpaths
    from . import scdg
except ImportError:
    pass

def get_tool(toolkey, func=False):
    """
    Parameters
    ----------
    func : bool, optional
        If True, return function, otherwise, return module.
    """
    tool = globals().get(toolkey)
    if func:
        return getattr(tool, toolkey)
    else:
        return tool

