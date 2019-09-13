import re
import sys
from datetime import datetime
from pathlib import Path

import black
import toml
from unidiff import PatchSet, PatchedFile

mode = black.FileMode(
    target_versions={black.TargetVersion.PY36},
    line_length=80,
    is_pyi=False,
    string_normalization=False,
)


def black_diff(src: Path) -> PatchedFile:
    src_contents = src.read_text()
    then = datetime.utcfromtimestamp(src.stat().st_mtime)
    now = datetime.utcnow()
    src_name = f"{src}\t{then} +0000"
    dst_name = f"{src}\t{now} +0000"
    try:
        dst_contents = black.format_file_contents(
            src_contents, fast=False, mode=mode
        )
        diff_str = black.diff(src_contents, dst_contents, src_name, dst_name)
        return PatchSet(diff_str.split('\n'))[0]
    except black.NothingChanged:
        return PatchedFile(src)


if __name__ == '__main__':
    thresh = int(sys.argv[1])
    ppt = toml.loads((Path() / 'pyproject.toml').read_text())
    exclude_re = re.compile(ppt['tool']['black']['exclude'], re.VERBOSE)

    excluded = [
        black_diff(src) for src in Path().glob('**/*.py')
        if exclude_re.match(str(Path('/') / src))
    ]
    for diff in sorted(excluded, key=lambda d: -(d.added + d.removed)):
        print(f'{diff.source_file}: +{diff.added} -{diff.removed}')
        if diff.added + diff.removed < thresh:
            print('File has < 10 changes to black formatting. Do it!')
            sys.exit(1)
    sys.exit(0)

    # def descend(root: Path, indent: int = 0):
    #     for file in root.iterdir():
    #         if file.is_dir():
    #             lines = list(descend(file, indent + 1))
    #             if not lines:
    #                 continue
    #             if len(lines) == 1:
    #                 yield f'{file.name}/{lines[0]}'
    #                 continue
    #             joined = f'\n{(indent+1) * "    "}|'.join(lines)
    #             yield (
    #                 f'{file.name}/(\n'
    #                 f'{indent * "    "}    {joined}\n'
    #                 f'{indent * "    "})'
    #             )
    #         elif file.suffix == '.py':
    #             diff = black_diff(file)
    #             if diff.added + diff.removed > thresh:
    #                 yield file.name
    #                 # yield f'{file.name}  # +{diff.added} -{diff.removed}'
    #
    # print("\n|".join(descend(Path())))
