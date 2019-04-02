## Contributing

Contributions to Scanpy are highly welcome!

### Before filing an issue

* Search the repository (also google) to see if someone has already reported the same issue. This allows contributors to spend less time responding to issues, and more time on adding features!
* Please provide a minimal complete verifiable example for any bug. If you're not sure what this means, check out [this blog post](http://matthewrocklin.com/blog/work/2018/02/28/minimal-bug-reports) by Matthew Rocklin or [this definition](https://stackoverflow.com/help/mcve) from StackOverflow.
* Let us know a bit about your environment. This can be as easy as pasting the results of `sc.logging.print_versions()`.

### Contributing code

* We stick to [PEP 8](https://www.python.org/dev/peps/pep-0008) and to this [editorconfig](https://github.com/theislab/scanpy/blob/master/.editorconfig) and *try* to stick to 80-character lines. In some cases, wider lines might improve readability, in most cases, not. Docstrings should always be 80 characters.
* We use the numpydoc style for writing docstrings. Either take a look at any Scanpy or Numpy function or [here](http://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_numpy.html). The `Params` abbreviation is a legit replacement for `Parameters`.
* Write tests for your functions! See [here](https://github.com/theislab/scanpy/tree/master/scanpy/tests) for examples.