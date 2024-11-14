#!/usr/bin/env python3
# /// script
# dependencies = [
#   "tomli; python_version < '3.11'",
#   "packaging",
# ]
# ///

from __future__ import annotations

import argparse
import sys
from collections import deque
from pathlib import Path
from typing import TYPE_CHECKING

if sys.version_info >= (3, 11):
    import tomllib
else:
    import tomli as tomllib

from packaging.requirements import Requirement
from packaging.version import Version

if TYPE_CHECKING:
    from collections.abc import Generator, Iterable


def min_dep(req: Requirement) -> Requirement:
    """
    Given a requirement, return the minimum version specifier.

    Example
    -------

    >>> min_dep(Requirement("numpy>=1.0"))
    <Requirement('numpy==1.0.*')>
    """
    req_name = req.name
    if req.extras:
        req_name = f"{req_name}[{','.join(req.extras)}]"

    filter_specs = [
        spec for spec in req.specifier if spec.operator in {"==", "~=", ">=", ">"}
    ]
    if not filter_specs:
        return Requirement(req_name)

    min_version = Version("0.0.0.a1")
    for spec in filter_specs:
        if spec.operator in {">", ">=", "~="}:
            min_version = max(min_version, Version(spec.version))
        elif spec.operator == "==":
            min_version = Version(spec.version)

    return Requirement(f"{req_name}=={min_version}.*")


def extract_min_deps(
    dependencies: Iterable[Requirement], *, pyproject
) -> Generator[Requirement, None, None]:
    dependencies = deque(dependencies)  # We'll be mutating this
    project_name = pyproject["project"]["name"]

    while len(dependencies) > 0:
        req = dependencies.pop()

        # If we are referring to other optional dependency lists, resolve them
        if req.name == project_name:
            assert req.extras, f"Project included itself as dependency, without specifying extras: {req}"
            for extra in req.extras:
                extra_deps = pyproject["project"]["optional-dependencies"][extra]
                dependencies += map(Requirement, extra_deps)
        else:
            yield min_dep(req)


def main():
    parser = argparse.ArgumentParser(
        prog="min-deps",
        description="""Parse a pyproject.toml file and output a list of minimum dependencies.

        Output is directly passable to `pip install`.""",
        usage="pip install `python min-deps.py pyproject.toml`",
    )
    parser.add_argument(
        "path", type=Path, help="pyproject.toml to parse minimum dependencies from"
    )
    parser.add_argument(
        "--extras", type=str, nargs="*", default=(), help="extras to install"
    )

    args = parser.parse_args()

    pyproject = tomllib.loads(args.path.read_text())

    project_name = pyproject["project"]["name"]
    deps = [
        *map(Requirement, pyproject["project"]["dependencies"]),
        *(Requirement(f"{project_name}[{extra}]") for extra in args.extras),
    ]

    min_deps = extract_min_deps(deps, pyproject=pyproject)

    print(" ".join(map(str, min_deps)))


if __name__ == "__main__":
    main()
