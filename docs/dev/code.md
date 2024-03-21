# Contributing code

## Development workflow

1. {ref}`Fork the Scanpy repository <forking-and-cloning>` to your own GitHub account
2. Create a {ref}`development environment <dev-environments>`
3. {ref}`Create a new branch <creating-a-branch>` for your PR
4. Add your feature or bugfix to the codebase
5. {ref}`Make sure all tests are passing <tests>`
6. {ref}`Build and visually check any changed documentation <building-the-docs>`
7. {ref}`Open a PR back to the main repository <open-a-pr>`

## Code style

New code should follow
[Black](https://black.readthedocs.io/en/stable/the_black_code_style.html)
and
[flake8](https://flake8.pycqa.org).
We ignore a couple of flake8 checks which are documented in the .flake8 file in the root of this repository.
To learn how to ignore checks per line please read
[flake8 violations](https://flake8.pycqa.org/en/latest/user/violations.html).
Additionally, we use Scanpy’s
[EditorConfig](https://github.com/scverse/scanpy/blob/main/.editorconfig),
so using an editor/IDE with support for both is helpful.
