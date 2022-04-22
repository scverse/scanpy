from contextlib import redirect_stdout
from datetime import datetime
from io import StringIO
import sys

import pytest

from scanpy import Verbosity, settings as s, logging as log
import scanpy as sc


@pytest.fixture
def logging_state():
    verbosity_orig = s.verbosity
    yield
    s.logfile = sys.stderr
    s.verbosity = verbosity_orig


def test_defaults():
    assert s.logpath is None


def test_formats(capsys, logging_state):
    s.logfile = sys.stderr
    s.verbosity = Verbosity.debug
    log.error('0')
    assert capsys.readouterr().err == 'ERROR: 0\n'
    log.warning('1')
    assert capsys.readouterr().err == 'WARNING: 1\n'
    log.info('2')
    assert capsys.readouterr().err == '2\n'
    log.hint('3')
    assert capsys.readouterr().err == '--> 3\n'
    log.debug('4')
    assert capsys.readouterr().err == '    4\n'


def test_deep(capsys, logging_state):
    s.logfile = sys.stderr
    s.verbosity = Verbosity.hint
    log.hint('0')
    assert capsys.readouterr().err == '--> 0\n'
    log.hint('1', deep='1!')
    assert capsys.readouterr().err == '--> 1\n'
    s.verbosity = Verbosity.debug
    log.hint('2')
    assert capsys.readouterr().err == '--> 2\n'
    log.hint('3', deep='3!')
    assert capsys.readouterr().err == '--> 3: 3!\n'


def test_logfile(tmp_path, logging_state):
    s.verbosity = Verbosity.hint

    io = StringIO()
    s.logfile = io
    assert s.logfile is io
    assert s.logpath is None
    log.error('test!')
    assert io.getvalue() == 'ERROR: test!\n'

    p = tmp_path / 'test.log'
    s.logpath = p
    assert s.logpath == p
    assert s.logfile.name == str(p)
    log.hint('test2')
    log.debug('invisible')
    assert s.logpath.read_text() == '--> test2\n'


def test_timing(monkeypatch, capsys, logging_state):
    s.logfile = sys.stderr
    counter = 0

    class IncTime:
        @staticmethod
        def now(tz):
            nonlocal counter
            counter += 1
            return datetime(2000, 1, 1, second=counter, microsecond=counter, tzinfo=tz)

    monkeypatch.setattr(log, 'datetime', IncTime)
    s.verbosity = Verbosity.debug

    log.hint('1')
    assert counter == 1 and capsys.readouterr().err == '--> 1\n'
    start = log.info('2')
    assert counter == 2 and capsys.readouterr().err == '2\n'
    log.hint('3')
    assert counter == 3 and capsys.readouterr().err == '--> 3\n'
    log.info('4', time=start)
    assert counter == 4 and capsys.readouterr().err == '4 (0:00:02)\n'
    log.info('5 {time_passed}', time=start)
    assert counter == 5 and capsys.readouterr().err == '5 0:00:03\n'


@pytest.mark.parametrize(
    "func",
    [
        sc.logging.print_header,
        sc.logging.print_versions,
        sc.logging.print_version_and_date,
    ],
)
def test_call_outputs(func):
    """
    Tests that these functions print to stdout and don't error.

    Checks that https://github.com/scverse/scanpy/issues/1437 is fixed.
    """
    output_io = StringIO()
    with redirect_stdout(output_io):
        func()
    output = output_io.getvalue()
    assert output != ""
