"""
Tools
"""

from . import diffmap
from . import diffrank
from . import dpt
from . import dbscan
from . import tsne
from . import sim
from . import spring

try:
    # development tools
    from . import paths
    from . import tgdyn
    from . import tgdyn_simple
except ImportError:
    pass

def get_tool(toolkey, func=False):
    """
    Parameters
    ----------
    func : bool, optional
        If True, return function, otherwise, return module.
    """
    tool_module = globals().get(toolkey)
    if func:
        return getattr(tool_module, toolkey)
    else:
        return tool_module

def get_plot_tool(toolkey, func=False):
    """
    Parameters
    ----------
    func : bool, optional
        If True, return function, otherwise, return module.
    """
    tool = globals().get(toolkey)
    if func:
        return getattr(tool, 'plot_' + toolkey)
    else:
        return tool

