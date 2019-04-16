import os
import sys
from argparse import ArgumentParser, Namespace, _SubParsersAction, ArgumentError
from collections import abc
from functools import lru_cache, partial
from pathlib import Path
from shutil import which
from subprocess import run, CompletedProcess
from typing import Optional, Generator, FrozenSet, Sequence, List, Tuple, Dict, Any, Mapping


class DelegatingSubparsersAction(_SubParsersAction):
    """Like a normal subcommand action, but uses a delegator for more choices"""
    def __init__(
        self,
        *args,
        _command: str,
        _runargs: Dict[str, Any],
        **kwargs,
    ):
        super().__init__(*args, **kwargs)
        self.command = _command
        self._name_parser_map = self.choices = CommandDelegator(_command, self, **_runargs)


class CommandDelegator(abc.MutableMapping):
    """Provide the ability to delegate, but don’t calculate the whole list until necessary"""
    def __init__(self, command: str, action: DelegatingSubparsersAction, **runargs):
        self.command = command
        self.action = action
        self.parser_map = {}
        self.runargs = runargs

    def __getitem__(self, k: str) -> ArgumentParser:
        try:
            return self.parser_map[k]
        except KeyError:
            if which(f'{self.command}-{k}'):
                return DelegatingParser(self, k)
            # Only here is the command list retrieved
            raise ArgumentError(self.action, f'No command “{k}”. Choose from {set(self)}')

    def __setitem__(self, k: str, v: ArgumentParser) -> None:
        self.parser_map[k] = v

    def __delitem__(self, k: str) -> None:
        del self.parser_map[k]

    # These methods retrieve the command list or help with doing it

    def __iter__(self) -> Generator[str, None, None]:
        yield from self.parser_map
        yield from self.commands

    def __len__(self) -> int:
        return len(self.parser_map) + len(self.commands)

    def __hash__(self) -> int:
        return hash(self.command)

    def __eq__(self, other: Mapping[str, ArgumentParser]):
        if isinstance(other, CommandDelegator):
            return all(
                getattr(self, attr) == getattr(other, attr)
                for attr in ['command', 'action', 'parser_map', 'runargs']
            )
        return self.parser_map == other

    @property
    @lru_cache()
    def commands(self) -> FrozenSet[str]:
        return frozenset(
            binary.name[len(self.command) + 1:]
            for bin_dir in os.environ['PATH'].split(os.pathsep)
            for binary in Path(bin_dir).glob(f'{self.command}-*')
            if os.access(binary, os.X_OK)
        )


class DelegatingParser(ArgumentParser):
    """Just sets parse_args().func to run the subcommand"""
    def __init__(self, cd: CommandDelegator, subcmd: str):
        super().__init__(f'{cd.command}-{subcmd}', add_help=False)
        self.cd = cd
        self.subcmd = subcmd

    def parse_known_args(
        self,
        args: Optional[Sequence[str]] = None,
        namespace: Optional[Namespace] = None,
    ) -> Tuple[Namespace, List[str]]:
        assert args is not None and namespace is None, 'Only use DelegatingParser as subparser'
        return Namespace(func=partial(run, [self.prog, *args], **self.cd.runargs)), []


def cmd_settings() -> None:
    from . import settings
    print(settings)


def main(argv: Optional[Sequence[str]] = None, *, check: bool = True, **runargs) -> Optional[CompletedProcess]:
    """Run a builtin scanpy command or a scanpy-* subcommand.

    Uses :func:`subcommand.run` for the latter: ``~run(['scanpy', *argv], **runargs)``
    """
    parser = ArgumentParser(description="There are a few packages providing commands. Try e.g. `pip install scanpy-scripts`!")
    parser.set_defaults(func=parser.print_help)
    
    subparsers: DelegatingSubparsersAction = parser.add_subparsers(
        action=DelegatingSubparsersAction,
        _command='scanpy',
        _runargs={**runargs, 'check': check},
    )
    
    parser_settings = subparsers.add_parser('settings')
    parser_settings.set_defaults(func=cmd_settings)
    
    args = parser.parse_args(argv)
    return args.func()


def console_main():
    """This serves as CLI entry point and will not show a Python traceback if a called command fails"""
    cmd = main(check=False)
    if cmd is not None:
        sys.exit(cmd.returncode)
