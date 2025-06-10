from __future__ import annotations

import os
from traceback import extract_stack
from types import FunctionType, MethodType


def documenting() -> bool:
    """Return whether this is being called from Sphinx."""
    if not os.environ.get("SPHINX_RUNNING"):
        return False
    for frame in extract_stack():
        # Let any sphinx ext get the docstring
        if frame.name in {
            "eval_config_file",  # Sphinx import
            "generate_autosummary_docs",  # Autosummary generator
            # "parse_generated_content",  # Autodoc parser
            "get_object_members",  # Class level of autodoc
            "import_object",  # Attr level of autodoc
        }:
            return True
    return False


class SingletonMeta(type):
    def __new__(mcls, cls_name: str, *args, **kwargs):
        cls = super().__new__(mcls, cls_name, *args, **kwargs)

        # We do something differently when we are imported by autosummary.
        if documenting():
            props = {}
            for name in dir(cls):
                if (attr := getattr(mcls, name, None)) is None:
                    continue
                if isinstance(attr, FunctionType | MethodType):
                    # Circumvent https://github.com/tox-dev/sphinx-autodoc-typehints/pull/157
                    setattr(cls, name, getattr(cls, name))
                if name not in cls.__dict__ and isinstance(attr, property):
                    # Allow autosummary to access the property, not the value
                    props[name] = getattr(mcls, name)

            def getattribute(_, name: str) -> object:
                """Return property or value depending on whether we are in autosummary.

                If an singleton instance property/method is accessed by autodoc/autosummary,
                return the property/method object, not the value/bound method.
                """
                if documenting() and name in props:
                    return props[name]
                return object.__getattribute__(cls, name)

            mcls.__getattribute__ = getattribute

        return cls

    def __dir__(cls) -> list[str]:
        # Deduplicate preserving order
        d = dict.fromkeys(super().__dir__()) | dict.fromkeys(dir(type(cls)))
        return [k for k in d if k != "mro"]
