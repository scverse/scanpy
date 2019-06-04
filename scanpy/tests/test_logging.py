import logging
import sys
from io import StringIO
from typing import Generator, List, Tuple

import pytest

from scanpy import Verbosity, settings as s, logging as l


@pytest.fixture
def logging_state():
    verbosity_orig = s.verbosity
    yield
    s.logfile = sys.stderr
    s.verbosity = verbosity_orig


def test_defaults():
    # This is only true when not in IPython,
    # But I hope we’ll never run tests from there
    assert s.logfile is sys.stderr
    assert s.logpath is None


def test_formats(logging_state):
    io = StringIO()

    def pop():
        v = io.getvalue()
        io.truncate(0)
        io.seek(0)
        return v

    s.verbosity = Verbosity.debug
    s.logfile = io
    l.error('0')
    assert pop() == 'ERROR: 0\n'
    l.warning('1')
    assert pop() == 'WARNING: 1\n'
    l.info('2')
    assert pop().endswith(' | 2\n')
    l.hint('3')
    assert pop() == '--> 3\n'
    l.debug('4')
    assert pop() == '    4\n'


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
