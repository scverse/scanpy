import os
import re
from pathlib import Path
from subprocess import PIPE
from typing import List

import pytest
from _pytest.capture import CaptureFixture
from _pytest.monkeypatch import MonkeyPatch

import scanpy
from scanpy.cli import main


HERE = Path(__file__).parent


@pytest.fixture
def set_path(monkeypatch: MonkeyPatch) -> None:
    monkeypatch.setenv('PATH', str(HERE / '_scripts'), prepend=os.pathsep)


def test_builtin_settings(capsys: CaptureFixture):
    main(['settings'])
    captured = capsys.readouterr()
    assert captured.out == f'{scanpy.settings}\n'


@pytest.mark.parametrize('args', [[], ['-h']])
def test_help_displayed(args: List[str], capsys: CaptureFixture):
    try:  # -h raises it, no args doesn’t. Maybe not ideal but meh.
        main(args)
    except SystemExit as se:
        assert se.code == 0
    captured = capsys.readouterr()
    assert captured.out.startswith('usage: ')


def test_help_output(set_path: type(None), capsys: CaptureFixture):
    with pytest.raises(SystemExit, match='^0$'):
        main(['-h'])
    captured = capsys.readouterr()
    assert re.search(
        r'^positional arguments:\n\s+\{settings,[\w,-]*testbin[\w,-]*\}$',
        captured.out,
        re.MULTILINE,
    )


def test_external(set_path: type(None)):
    # We need to capture the output manually, since subprocesses don’t write to sys.stderr
    cmdline = ['testbin', '-t', '--testarg', 'testpos']
    cmd = main(cmdline, stdout=PIPE, encoding='utf-8', check=True)
    assert cmd.stdout == 'test -t --testarg testpos\n'


def test_error_wrong_command(capsys: CaptureFixture):
    with pytest.raises(SystemExit, match='^2$'):
        main(['idonotexist--'])
    captured = capsys.readouterr()
    assert 'No command “idonotexist--”. Choose from' in captured.err
