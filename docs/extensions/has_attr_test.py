from jinja2.defaults import DEFAULT_NAMESPACE
from jinja2.utils import import_string
from sphinx.application import Sphinx


def has_attr(obj_path: str, attr: str) -> bool:
    # https://jinja.palletsprojects.com/en/3.0.x/api/#custom-tests
    obj = import_string(obj_path)
    return hasattr(obj, attr)


def setup(app: Sphinx):
    DEFAULT_NAMESPACE["has_attr"] = has_attr
