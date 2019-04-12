from argparse import ArgumentParser
from collections import abc
from typing import Mapping, Generator

parser = ArgumentParser(description="There are a few packages providing commands. Try e.g. `pip install scanpy-scripts`!")

subparsers = parser.add_subparsers()

subparsers.add_parser('settings')


class Delegator(abc.Mapping):
    def __init__(self, parser_map: Mapping[str, ArgumentParser]):
        self.parser_map = parser_map

    def __iter__(self) -> Generator[str, None, None]:
        yield from self.parser_map
        yield from self._find_commands()

    def __len__(self) -> int:
        raise NotImplementedError

    def __getitem__(self, k: str) -> ArgumentParser:
        raise NotImplementedError

    def _find_commands(self) -> Generator[str, None, None]:
        raise NotImplementedError


# has to happen last
subparsers._name_parser_map = Delegator(subparsers._name_parser_map)

args = parser.parse_args()

# TODO
