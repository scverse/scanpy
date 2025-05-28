#!/usr/bin/env python3
# /// script
# dependencies = [
#   "tomli; python_version < '3.11'",
#   "packaging",
# ]
# ///
"""Parse a pyproject.toml file and output a list of minimum dependency versions."""

from __future__ import annotations

import argparse
import sys
from collections import deque
from contextlib import ExitStack
from functools import cached_property
from pathlib import Path
from typing import TYPE_CHECKING

if sys.version_info >= (3, 11):
    import tomllib
else:
    import tomli as tomllib

from packaging.requirements import Requirement
from packaging.version import Version

if TYPE_CHECKING:
    from collections.abc import Generator, Iterable, Sequence
    from collections.abc import Set as AbstractSet
    from typing import Any, Self


def min_dep(req: Requirement) -> Requirement:
    """Given a requirement, return the minimum version specifier.

    Example
    -------
    >>> min_dep(Requirement("numpy>=1.0"))
    <Requirement('numpy==1.0')>
    >>> min_dep(Requirement("numpy<3.0"))
    <Requirement('numpy<3.0')>
    """
    req_name = req.name
    if req.extras:
        req_name = f"{req_name}[{','.join(req.extras)}]"

    filter_specs = [
        spec for spec in req.specifier if spec.operator in {"==", "~=", ">=", ">"}
    ]
    if not filter_specs:
        # TODO: handle markers
        return Requirement(f"{req_name}{req.specifier}")
    min_version = Version("0.0.0.a1")
    for spec in filter_specs:
        if spec.operator in {">", ">=", "~="}:
            min_version = max(min_version, Version(spec.version))
        elif spec.operator == "==":
            min_version = Version(spec.version)

    return Requirement(f"{req_name}=={min_version}")


def extract_min_deps(
    dependencies: Iterable[Requirement], *, pyproject
) -> Generator[Requirement, None, None]:
    """Extract minimum dependency versions from a list of requirements."""
    dependencies = deque(dependencies)  # We'll be mutating this
    project_name = pyproject["project"]["name"]

    deps = {}
    while len(dependencies) > 0:
        req = dependencies.pop()

        # If we are referring to other optional dependency lists, resolve them
        if req.name == project_name:
            assert req.extras, (
                f"Project included itself as dependency, without specifying extras: {req}"
            )
            for extra in req.extras:
                extra_deps = pyproject["project"]["optional-dependencies"][extra]
                dependencies += map(Requirement, extra_deps)
        else:
            if req.name in deps:
                req.specifier &= deps[req.name].specifier
                req.extras |= deps[req.name].extras
            deps[req.name] = min_dep(req)
    yield from deps.values()


class Args(argparse.Namespace):
    """Parse a pyproject.toml file and output a list of minimum dependencies.

    Output is optimized for `[uv] pip install` (see `-o`/`--output` for details).
    """

    _path: Path
    output: Path | None
    _extras: list[str]
    _all_extras: bool

    @classmethod
    def parse(cls, argv: Sequence[str] | None = None) -> Self:
        """Parse CLI arguments."""
        return cls.parser().parse_args(argv, cls())

    @classmethod
    def parser(cls) -> argparse.ArgumentParser:
        """Construct a CLI argument parser."""
        parser = argparse.ArgumentParser(
            prog="min-deps",
            description=cls.__doc__,
            usage="pip install `python min-deps.py pyproject.toml`",
            allow_abbrev=False,
        )
        parser.add_argument(
            "_path",
            metavar="pyproject.toml",
            type=Path,
            help="Path to pyproject.toml to parse minimum dependencies from",
        )
        parser.add_argument(
            "--extras",
            dest="_extras",
            metavar="EXTRA",
            type=str,
            nargs="*",
            default=(),
            help="extras to install",
        )
        parser.add_argument(
            "--all-extras",
            dest="_all_extras",
            action="store_true",
            help="get all extras",
        )
        parser.add_argument(
            *("--output", "-o"),
            metavar="FILE",
            type=Path,
            default=None,
            help=(
                "output file (default: stdout). "
                "Without this option, output is space-separated for direct passing to `pip install`. "
                "With this option, output written to a file newline-separated file usable as `requirements.txt` or `constraints.txt`."
            ),
        )
        return parser

    @cached_property
    def pyproject(self) -> dict[str, Any]:
        """Return the parsed `pyproject.toml`."""
        return tomllib.loads(self._path.read_text())

    @cached_property
    def extras(self) -> AbstractSet[str]:
        """Return the extras to install."""
        if self._extras:
            if self._all_extras:
                sys.exit("Cannot specify both --extras and --all-extras")
            return dict.fromkeys(self._extras).keys()
        if not self._all_extras:
            return set()
        return self.pyproject["project"]["optional-dependencies"].keys()


def main(argv: Sequence[str] | None = None) -> None:
    """Run main entry point."""
    args = Args.parse(argv)

    project_name = args.pyproject["project"]["name"]
    deps = [
        *map(Requirement, args.pyproject["project"]["dependencies"]),
        *(Requirement(f"{project_name}[{extra}]") for extra in args.extras),
    ]

    min_deps = extract_min_deps(deps, pyproject=args.pyproject)

    sep = "\n" if args.output else " "
    with ExitStack() as stack:
        f = stack.enter_context(args.output.open("w")) if args.output else sys.stdout
        print(sep.join(map(str, min_deps)), file=f)


if __name__ == "__main__":
    main()
