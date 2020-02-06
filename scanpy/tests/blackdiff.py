import re
import sys
from collections import Counter
from difflib import Differ
from pathlib import Path
from typing import Tuple

import black
import toml

mode = black.FileMode(
    target_versions={black.TargetVersion.PY36},
    line_length=80,
    is_pyi=False,
    string_normalization=False,
)


def black_diff(src: Path, differ=Differ()) -> Tuple[int, int]:
    src_contents = src.read_text()
    dst_contents = black.format_str(src_contents, mode=mode)
    if src_contents == dst_contents:
        return 0, 0

    counts = Counter(
        line[0]
        for line in differ.compare(src_contents.splitlines(), dst_contents.splitlines())
    )
    return counts['+'], counts['-']


if __name__ == '__main__':
    thresh = int(sys.argv[1])
    ppt = toml.loads((Path() / 'pyproject.toml').read_text())
    exclude_re = re.compile(ppt['tool']['black']['exclude'], re.VERBOSE)

    excluded = [
        (src, black_diff(src))
        for src in Path().glob('**/*.py')
        if exclude_re.match(str(Path('/') / src))
        and not src.parts[0] == 'build'
        and not src.parts[:2] == ('scanpy', 'api')
    ]
    for file, (added, removed) in sorted(excluded, key=lambda sd: -sum(sd[1])):
        print(f'{file}: +{added} -{removed}')
        if added == removed == 0:
            print(
                'This above file has 0 changes to black formatting, good job! '
                'Please remove it from “tool.black.exclude" in pyproject.toml'
            )
            sys.exit(1)
        elif added + removed < thresh:
            print(
                f'This above file has < {thresh} changes to black formatting. '
                'Please remove it from “tool.black.exclude" in pyproject.toml'
                'and afterwards black format it.'
            )
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
