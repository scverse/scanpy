from __future__ import annotations

import sys
from contextlib import redirect_stdout
from datetime import datetime
from io import StringIO
from logging import StreamHandler
from typing import TYPE_CHECKING

import pytest

import scanpy as sc
from scanpy import Verbosity
from scanpy import logging as log
from scanpy import settings as s

if TYPE_CHECKING:
    from collections.abc import Mapping
    from pathlib import Path


def test_defaults(
    caplog: pytest.LogCaptureFixture, original_settings: Mapping[str, object]
) -> None:
    assert s.logpath is original_settings["_logpath"] is None
    assert s.logfile is original_settings["_logfile"] is sys.stderr
    # we override s.verbosity, so we only check the default here:
    assert original_settings["_verbosity"] is Verbosity.warning

    # check logging handler file and level
    [handler] = (h for h in s._root_logger.handlers if h is not caplog.handler)
    assert isinstance(handler, StreamHandler)
    assert handler.stream is s.logfile
    assert s._root_logger.level == s.verbosity.level


def test_records(caplog: pytest.LogCaptureFixture) -> None:
    s.verbosity = Verbosity.debug
    log.error("0")
    log.warning("1")
    log.info("2")
    log.hint("3")
    log.debug("4")
    assert caplog.record_tuples == [
        ("root", 40, "0"),
        ("root", 30, "1"),
        ("root", 20, "2"),
        ("root", 15, "3"),
        ("root", 10, "4"),
    ]


def test_formats(capsys: pytest.CaptureFixture):
    s.logfile = sys.stderr
    s.verbosity = Verbosity.debug
    log.error("0")
    assert capsys.readouterr().err == "ERROR: 0\n"
    log.warning("1")
    assert capsys.readouterr().err == "WARNING: 1\n"
    log.info("2")
    assert capsys.readouterr().err == "2\n"
    log.hint("3")
    assert capsys.readouterr().err == "--> 3\n"
    log.debug("4")
    assert capsys.readouterr().err == "    4\n"


def test_deep(capsys: pytest.CaptureFixture):
    s.logfile = sys.stderr
    s.verbosity = Verbosity.hint
    log.hint("0")
    assert capsys.readouterr().err == "--> 0\n"
    log.hint("1", deep="1!")
    assert capsys.readouterr().err == "--> 1\n"
    s.verbosity = Verbosity.debug
    log.hint("2")
    assert capsys.readouterr().err == "--> 2\n"
    log.hint("3", deep="3!")
    assert capsys.readouterr().err == "--> 3: 3!\n"


def test_logfile(tmp_path: Path, caplog: pytest.LogCaptureFixture):
    s.verbosity = Verbosity.hint

    io = StringIO()
    s.logfile = io
    assert s.logfile is io
    assert s.logpath is None
    log.error("test!")
    assert io.getvalue() == "ERROR: test!\n"

    # setting a logfile removes all handlers
    assert not caplog.records

    p = tmp_path / "test.log"
    s.logpath = p
    assert s.logpath == p
    assert s.logfile.name == str(p)
    log.hint("test2")
    log.debug("invisible")
    assert s.logpath.read_text() == "--> test2\n"

    # setting a logfile removes all handlers
    assert not caplog.records


def test_timing(monkeypatch, capsys: pytest.CaptureFixture):
    counter = 0

    class IncTime:
        @staticmethod
        def now(tz):
            nonlocal counter
            counter += 1
            return datetime(2000, 1, 1, second=counter, microsecond=counter, tzinfo=tz)

    monkeypatch.setattr(log, "datetime", IncTime)
    s.logfile = sys.stderr
    s.verbosity = Verbosity.debug

    log.hint("1")
    assert counter == 1
    assert capsys.readouterr().err == "--> 1\n"

    start = log.info("2")
    assert counter == 2
    assert capsys.readouterr().err == "2\n"

    log.hint("3")
    assert counter == 3
    assert capsys.readouterr().err == "--> 3\n"

    log.info("4", time=start)
    assert counter == 4
    assert capsys.readouterr().err == "4 (0:00:02)\n"

    log.info("5 {time_passed}", time=start)
    assert counter == 5
    assert capsys.readouterr().err == "5 0:00:03\n"


@pytest.mark.parametrize(
    "func",
    [
        sc.logging.print_header,
        sc.logging.print_versions,
        sc.logging.print_version_and_date,
    ],
)
def test_call_outputs(func):
    """Tests that these functions print to stdout and don't error.

    Checks that https://github.com/scverse/scanpy/issues/1437 is fixed.
    """
    output_io = StringIO()
    with redirect_stdout(output_io):
        out = func()
        if out is not None:
            print(out)
    output = output_io.getvalue()
    assert output != ""
