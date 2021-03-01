"""Temporary setuptools bridge
Don't use this except if you have a deadline or you encounter a bug.
"""
import re
import sys
from warnings import warn

import setuptools
from pathlib import Path

from flit_core import common, config
from setuptools_scm.integration import find_files


field_map = dict(
    description="summary",
    long_description="description",
    long_description_content_type="description_content_type",
    python_requires="requires_python",
    url="home_page",
    **{
        n: n
        for n in ["name", "version", "author", "author_email", "license", "classifiers"]
    },
)


def setup_args(config_path=Path("pyproject.toml")):
    cfg = config.read_flit_config(config_path)
    module = common.Module(cfg.module, config_path.parent)
    metadata = common.make_metadata(module, cfg)
    kwargs = {}
    for st_field, metadata_field in field_map.items():
        val = getattr(metadata, metadata_field, None)
        if val is not None:
            kwargs[st_field] = val
        else:
            msg = f"{metadata_field} not found in {dir(metadata)}"
            assert metadata_field in {"license"}, msg
    kwargs["packages"] = setuptools.find_packages(include=[metadata.name + "*"])
    if metadata.requires_dist:
        kwargs["install_requires"] = [
            req for req in metadata.requires_dist if "extra ==" not in req
        ]
    if cfg.reqs_by_extra:
        kwargs["extras_require"] = cfg.reqs_by_extra
    scripts = cfg.entrypoints.get("console_scripts")
    if scripts is not None:
        kwargs["entry_points"] = dict(
            console_scipts=[" = ".join(ep) for ep in scripts.items()]
        )
    kwargs["include_package_data"] = True
    kwargs["package_data"] = {
        module.name: [
            re.escape(f[len(module.name) + 1 :]) for f in find_files(module.path)
        ]
    }
    return kwargs


if __name__ == "__main__":
    if "develop" in sys.argv:
        msg = (
            "Please use `flit install -s` or `flit install --pth-file` "
            "instead of `pip install -e`/`python setup.py develop` "
            "once https://github.com/pypa/pip/issues/9670 is fixed."
        )
    elif "install" in sys.argv:
        msg = 'Please use `pip install "$d"` instead of `python "$d/setup.py" install`'
    else:
        msg = "Please use `pip ...` or `flit ...` instead of `python setup.py ...`"
    warn(msg, FutureWarning)
    setuptools.setup(**setup_args())
