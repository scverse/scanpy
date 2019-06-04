import sys
from io import StringIO

from scanpy import Verbosity, settings as s, logging as l


def test_defaults():
    pass


def test_logfile(tmp_path):
    s.verbosity = Verbosity.hint
    assert s.logfile is sys.stderr
    assert s.logpath is None

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
