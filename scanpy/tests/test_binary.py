import os
from pathlib import Path
from subprocess import PIPE

import pytest
from _pytest.capture import CaptureFixture
from _pytest.monkeypatch import MonkeyPatch

import scanpy
from scanpy.cli import main


HERE = Path(__file__).parent


@pytest.fixture
def set_path(monkeypatch: MonkeyPatch):
    monkeypatch.setenv('PATH', str(HERE / '_scripts'), prepend=os.pathsep)


def test_help_no_command(capsys: CaptureFixture):
    main([])
    captured = capsys.readouterr()
    assert captured.out.startswith('usage: ')


def test_builtin_settings(capsys: CaptureFixture):
    main(['settings'])
    captured = capsys.readouterr()
    assert captured.out == f'{scanpy.settings}\n'


def test_external(set_path):
    # We need to capture the output manually, since subprocesses don’t write to sys.stderr
    cmd = main(['testbin', '-t', '--testarg', 'testpos'], stdout=PIPE, text=True, check=True)
    assert cmd.stdout == 'test -t --testarg testpos\n'


def test_error_wrong_command(capsys: CaptureFixture):
    with pytest.raises(SystemExit, match=r'^2$'):
        main(['idonotexist--'])
    captured = capsys.readouterr()
    assert 'No command “idonotexist--”. Choose from' in captured.err
