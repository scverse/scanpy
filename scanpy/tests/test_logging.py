import sys
from datetime import datetime
from io import StringIO

import pytest

from scanpy import Verbosity, settings as s, logging as l


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
    l.error('0')
    assert capsys.readouterr().err == 'ERROR: 0\n'
    l.warning('1')
    assert capsys.readouterr().err == 'WARNING: 1\n'
    l.info('2')
    assert capsys.readouterr().err == '2\n'
    l.hint('3')
    assert capsys.readouterr().err == '--> 3\n'
    l.debug('4')
    assert capsys.readouterr().err == '    4\n'


def test_deep(capsys, logging_state):
    s.logfile = sys.stderr
    s.verbosity = Verbosity.hint
    l.hint('0')
    assert capsys.readouterr().err == '--> 0\n'
    l.hint('1', deep='1!')
    assert capsys.readouterr().err == '--> 1\n'
    s.verbosity = Verbosity.debug
    l.hint('2')
    assert capsys.readouterr().err == '--> 2\n'
    l.hint('3', deep='3!')
    assert capsys.readouterr().err == '--> 3: 3!\n'


def test_logfile(tmp_path, logging_state):
    s.verbosity = Verbosity.hint

    io = StringIO()
    s.logfile = io
    assert s.logfile is io
    assert s.logpath is None
    l.error('test!')
    assert io.getvalue() == 'ERROR: test!\n'

    p = tmp_path / 'test.log'
    s.logpath = p
    assert s.logpath == p
    assert s.logfile.name == str(p)
    l.hint('test2')
    l.debug('invisible')
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

    monkeypatch.setattr(l, 'datetime', IncTime)
    s.verbosity = Verbosity.debug

    l.hint('1')
    assert counter == 1 and capsys.readouterr().err == '--> 1\n'
    start = l.info('2')
    assert counter == 2 and capsys.readouterr().err == '2\n'
    l.hint('3')
    assert counter == 3 and capsys.readouterr().err == '--> 3\n'
    l.info('4', time=start)
    assert counter == 4 and capsys.readouterr().err == '4 (0:00:02)\n'
    l.info('5 {time_passed}', time=start)
    assert counter == 5 and capsys.readouterr().err == '5 0:00:03\n'
