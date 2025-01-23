# Contributing code

## Development workflow

1. {ref}`Fork the Scanpy repository <forking-and-cloning>` to your own GitHub account
2. Create a {ref}`development environment <dev-environments>`
3. {ref}`Create a new branch <creating-a-branch>` for your PR
4. Add your feature or bugfix to the codebase
5. {ref}`Make sure all tests are passing <tests>`
6. {ref}`Build and visually check any changed documentation <building-the-docs>`
7. {ref}`Open a PR back to the main repository <open-a-pr>`
8. {ref}`Add a release note to your PR <adding-to-the-docs>`

## Code style

Code contributions will be formatted and style checked using [Ruff][].
Ignored checks are configured in the `tool.ruff.lint` section of {file}`pyproject.toml`.
To learn how to ignore checks per line please read about [ignoring errors][].
Additionally, we use Scanpy’s [EditorConfig][],
so using an editor/IDE with support for both is helpful.

[Ruff]: https://docs.astral.sh/ruff/
[ignoring errors]: https://docs.astral.sh/ruff/tutorial/#ignoring-errors
[EditorConfig]: https://github.com/scverse/scanpy/blob/main/.editorconfig
