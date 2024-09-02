from __future__ import annotations

import os
import re
from contextlib import nullcontext
from pathlib import Path
from subprocess import PIPE
from typing import TYPE_CHECKING

import pytest

import scanpy
from scanpy.cli import main

if TYPE_CHECKING:
    from _pytest.capture import CaptureFixture
    from _pytest.monkeypatch import MonkeyPatch

HERE = Path(__file__).parent


@pytest.fixture
def _set_path(monkeypatch: MonkeyPatch) -> None:
    monkeypatch.setenv("PATH", str(HERE / "_scripts"), prepend=os.pathsep)


def test_builtin_settings(capsys: CaptureFixture):
    main(["settings"])
    captured = capsys.readouterr()
    assert captured.out == f"{scanpy.settings}\n"


@pytest.mark.parametrize("args", [[], ["-h"]])
def test_help_displayed(args: list[str], capsys: CaptureFixture):
    # -h raises it, no args doesn’t. Maybe not ideal but meh.
    ctx = pytest.raises(SystemExit) if args else nullcontext()
    with ctx as se:
        main(args)
    if se is not None:
        assert se.value.code == 0
    captured = capsys.readouterr()
    assert captured.out.startswith("usage: ")


@pytest.mark.usefixtures("_set_path")
def test_help_output(capsys: CaptureFixture):
    with pytest.raises(SystemExit, match="^0$"):
        main(["-h"])
    captured = capsys.readouterr()
    assert re.search(
        r"^positional arguments:\n\s+\{settings,[\w,-]*testbin[\w,-]*\}$",
        captured.out,
        re.MULTILINE,
    )


@pytest.mark.usefixtures("_set_path")
def test_external():
    # We need to capture the output manually, since subprocesses don’t write to sys.stderr
    cmdline = ["testbin", "-t", "--testarg", "testpos"]
    cmd = main(cmdline, stdout=PIPE, encoding="utf-8", check=True)
    assert cmd.stdout == "test -t --testarg testpos\n"


def test_error_wrong_command(capsys: CaptureFixture):
    with pytest.raises(SystemExit, match="^2$"):
        main(["idonotexist--"])
    captured = capsys.readouterr()
    assert "invalid choice: 'idonotexist--' (choose from" in captured.err
