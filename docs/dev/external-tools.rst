External tools
==============

For external tools we’d prefer the code mostly lives in an external package.
For internal, it’ll be a more in-depth review as we'd be more responsible for maintenance and would essentially be endorsing the method.

Docs
----

To make sure your tool gets an entry in the documentation, you'll have to add it.
Make sure there is a reference to your new function in `docs/external.rst`.

Testing
-------

Tests for tools in the external module are less thorough than tools in the main API.

While tests for the main API are there to assure correctness, we assume external tools have their own test suites in their package.
Instead, tests here are to let us know:

* If changes in scanpy break the external tool
* If changes in dependencies would break external tools

If the tool is not installed, the test should not make the suite fail.
You can do this by starting the tests with: :func:`~pytest.importorskip`.

To make sure the test runs on CI, you'll need to make sure it gets installed there.
This requires:

* An entry to the `extras_requires` field in `setup.py`
* Modifying the installation command in the CI configs (e.g. `.azure-pipelines.yml`).
