"""Definition for scanpy’s CLI entry point to be used programmatically."""

from __future__ import annotations

import os
import sys
from argparse import ArgumentParser, Namespace, _SubParsersAction
from collections.abc import MutableMapping
from functools import cached_property, partial
from pathlib import Path
from shutil import which
from subprocess import run
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Iterator, Mapping, Sequence
    from subprocess import CompletedProcess
    from typing import Any


class _DelegatingSubparsersAction(_SubParsersAction):
    """Like a normal subcommand action, but uses a delegator for more choices."""

    def __init__(self, *args, _command: str, _runargs: dict[str, Any], **kwargs):
        super().__init__(*args, **kwargs)
        self.command = _command
        self._name_parser_map = self.choices = _CommandDelegator(
            _command, self, **_runargs
        )


class _CommandDelegator(MutableMapping):
    """Provide the ability to delegate, but don’t calculate the whole list until necessary."""

    def __init__(self, command: str, action: _DelegatingSubparsersAction, **runargs):
        self.command = command
        self.action = action
        self.parser_map = {}
        self.runargs = runargs

    def __contains__(self, k: str) -> bool:
        if k in self.parser_map:
            return True
        try:
            self[k]
        except KeyError:
            return False
        return True

    def __getitem__(self, k: str) -> ArgumentParser:
        try:
            return self.parser_map[k]
        except KeyError:
            if which(f"{self.command}-{k}"):
                return _DelegatingParser(self, k)
            raise

    def __setitem__(self, k: str, v: ArgumentParser) -> None:
        self.parser_map[k] = v

    def __delitem__(self, k: str) -> None:
        del self.parser_map[k]

    # These methods retrieve the command list or help with doing it

    def __iter__(self) -> Iterator[str]:
        yield from self.parser_map
        yield from self.commands

    def __len__(self) -> int:
        return len(self.parser_map) + len(self.commands)

    def __hash__(self) -> int:
        return hash(self.command)

    def __eq__(self, other: Mapping[str, ArgumentParser]):
        if isinstance(other, _CommandDelegator):
            return all(
                getattr(self, attr) == getattr(other, attr)
                for attr in ["command", "action", "parser_map", "runargs"]
            )
        return self.parser_map == other

    @cached_property
    def commands(self) -> frozenset[str]:
        return frozenset(
            binary.name[len(self.command) + 1 :]
            for bin_dir in os.environ["PATH"].split(os.pathsep)
            for binary in Path(bin_dir).glob(f"{self.command}-*")
            if os.access(binary, os.X_OK)
        )


class _DelegatingParser(ArgumentParser):
    """Just sets parse_args().func to run the subcommand."""

    def __init__(self, cd: _CommandDelegator, subcmd: str):
        super().__init__(f"{cd.command}-{subcmd}", add_help=False)
        self.cd = cd
        self.subcmd = subcmd

    def parse_known_args(
        self,
        args: Sequence[str] | None = None,
        namespace: Namespace | None = None,
    ) -> tuple[Namespace, list[str]]:
        msg = "Only use DelegatingParser as subparser"
        assert args is not None, msg
        assert namespace is None, msg
        return Namespace(func=partial(run, [self.prog, *args], **self.cd.runargs)), []


def _cmd_settings() -> None:
    from ._settings import settings

    print(settings)


def main(
    argv: Sequence[str] | None = None, *, check: bool = True, **runargs
) -> CompletedProcess | None:
    """Run a builtin scanpy command or a scanpy-* subcommand.

    Uses :func:`subcommand.run` for the latter:
    `~run(['scanpy', *argv], **runargs)`
    """
    parser = ArgumentParser(
        description=(
            "There are a few packages providing commands. "
            "Try e.g. `pip install scanpy-scripts`!"
        )
    )
    parser.set_defaults(func=parser.print_help)

    subparsers: _DelegatingSubparsersAction = parser.add_subparsers(
        action=_DelegatingSubparsersAction,
        _command="scanpy",
        _runargs={**runargs, "check": check},
    )

    parser_settings = subparsers.add_parser("settings")
    parser_settings.set_defaults(func=_cmd_settings)

    args = parser.parse_args(argv)
    return args.func()


def console_main():
    """Serve as CLI entry point and don’t show a Python traceback if a called command fails."""
    cmd = main(check=False)
    if cmd is not None:
        sys.exit(cmd.returncode)
