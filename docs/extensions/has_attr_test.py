"""Extension adding a jinja2 filter that tests if an object has an attribute."""

from __future__ import annotations

from inspect import get_annotations
from typing import TYPE_CHECKING

from jinja2.defaults import DEFAULT_NAMESPACE
from jinja2.utils import import_string

if TYPE_CHECKING:
    from sphinx.application import Sphinx


def has_member(obj_path: str, attr: str) -> bool:
    """Test if an object has an attribute."""
    # https://jinja.palletsprojects.com/en/3.0.x/api/#custom-tests
    obj = import_string(obj_path)
    return hasattr(obj, attr) or attr in get_annotations(obj)


def setup(app: Sphinx):
    """App setup hook."""
    DEFAULT_NAMESPACE["has_member"] = has_member
