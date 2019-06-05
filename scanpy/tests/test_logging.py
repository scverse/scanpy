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


class Recorder:
    def __init__(self):
        self.records = []

    def write(self, s: str):
        self.records.append(s)

    def pop(self):
        return self.records.pop()


def test_defaults():
    # This is only true when not in IPython,
    # But I hope we’ll never run tests from there
    assert s.logfile is sys.stderr
    assert s.logpath is None


def test_formats(logging_state):
    s.verbosity = Verbosity.debug
    rec = Recorder()
    s.logfile = rec
    l.error('0')
    assert rec.pop() == 'ERROR: 0\n'
    l.warning('1')
    assert rec.pop() == 'WARNING: 1\n'
    l.info('2')
    assert rec.pop().endswith(' | 2\n')
    l.hint('3')
    assert rec.pop() == '--> 3\n'
    l.debug('4')
    assert rec.pop() == '    4\n'


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


def test_timing(monkeypatch, logging_state):
    counter = 0

    class IncTime:
        @staticmethod
        def now(tz):
            nonlocal counter
            counter += 1
            return datetime(2000, 1, 1, second=counter, tzinfo=tz)

    monkeypatch.setattr(l, 'datetime', IncTime)
    s.verbosity = Verbosity.debug
    rec = Recorder()
    s.logfile = rec

    l.hint('1')
    assert counter == 1 and rec.pop() == '--> 1\n'
    start = l.info('2')
    assert counter == 2 and rec.pop().endswith(' | 2\n')
    l.hint('3')
    assert counter == 3 and rec.pop() == '--> 3\n'
    l.info('4', time=start)
    assert counter == 4 and rec.pop().endswith(' | 4 (0:00:02)\n')
