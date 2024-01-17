#!python3
from __future__ import annotations

import argparse
import sys
from pathlib import Path

if sys.version_info >= (3, 11):
    import tomllib
else:
    import tomli as tomllib
from packaging.requirements import Requirement
from packaging.version import Version


def min_dep(req: Requirement) -> str:
    """
    Given a requirement, return the minimum version specifier.

    Example
    -------

    >>> min_dep(Requirement("numpy>=1.0"))
    "numpy==1.0"
    """
    req_name = req.name
    if req.extras:
        req_name = f"{req_name}[{','.join(req.extras)}]"

    # TODO: Should this be allowed?
    if not req.specifier:
        return req_name

    min_version = Version("0.0.0.a1")
    for spec in req.specifier:
        if spec.operator in [">", ">=", "~-"]:
            min_version = max(min_version, Version(spec.version))
        elif spec.operator == "==":
            min_version = Version(spec.version)

    # TODO: should this return `~=` or `==`?
    return f"{req_name}=={min_version}.*"


def extract_min_deps(
        dependencies: list[str], 
        *, 
        pyproject
    ) -> list[str]:
    dependencies = dependencies.copy()  # We'll be mutating this
    requirements: list[Requirement] = []
    project_name = pyproject["project"]["name"]

    while len(dependencies) > 0:
        req = Requirement(dependencies.pop())

        # If we are reffering to other optional dependency lists, resolve them
        if req.name == project_name:
            assert req.extras, f"Project included itself as dependency, without specifying extras: {req}"
            for extra in req.extras:
                dependencies.extend(pyproject["project"]["optional-dependencies"][extra])
        else:
            requirements.append(min_dep(req))

    return requirements


def main():
    # TODO: Allow optional dependencies
    parser = argparse.ArgumentParser(
        prog="min-deps",
        description="""Parse a pyproject.toml file and output a list of minimum dependencies.

        Output is directly passable to `pip install`.""",
        usage="pip install `python min-deps.py pyproject.toml`",
    )
    parser.add_argument(
        "path", type=Path, help="pyproject.toml to parse minimum dependencies from"
    )
    parser.add_argument("--extras", type=str, nargs="*", help="extras to install")

    args = parser.parse_args()

    pyproject = tomllib.loads(args.path.read_text())

    project_name = pyproject["project"]["name"]
    deps = pyproject["project"]["dependencies"]

    for extra in args.extras:
        deps.append(f"{project_name}[{extra}]")

    min_deps = extract_min_deps(deps, pyproject=pyproject)

    print(" ".join(min_deps))


if __name__ == "__main__":
    main()
