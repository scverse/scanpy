import sys

from .cli import main


cmd = main(check=False)

if cmd is not None:
    sys.exit(cmd.returncode)
