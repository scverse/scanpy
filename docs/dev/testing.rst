Tests
=====

Possibly the most important part of contributing to any open source package is the test suite.
Implementations may change, but the only way we can know the code is working before making a release is the test suite.

Running the tests
-----------------

We use `pytest <https://docs.pytest.org/en/stable/>`__ to test scanpy.
To run the tests first make sure you have the required dependencies (``pip install -e "scanpy[tests]"``), then run ``pytest`` from the root of the repository.

It can take a while to run the whole test suite. There are a few ways to cut down on this while working on a PR:

1. Only run a subset of the tests. This can be done with the `-k` argument from pytest (e.g. ``pytest -k test_plotting.py`` or ``pytest -k "test_umap*"``
2. Run tests in parallel. If you install the pytest extension `pytest-xdist <https://github.com/pytest-dev/pytest-xdist>`__ you can run tests in parallel with the ``--numprocesses`` argument to pytest (e.g. ``pytest -n 8``).

A lot of warnings can be thrown while running the test suite. It's often easier to parse the output with these excluded (

Writing tests
-------------

You can refer to the `existing test suite <https://github.com/theislab/scanpy/tree/master/scanpy/tests>`__ for examples.
If you haven't written tests before, Software Carpentry has an `in-depth guide <http://katyhuff.github.io/python-testing/>`__ on the topic.

Some tips on writing tests:

We highly recommend using `Test Driven Development <https://en.wikipedia.org/wiki/Test-driven_development>`__ when contributing code.
This process boils down to writing the tests before you start the implementation.


Performance
~~~~~~~~~~~

It's more important that you're accurately testing the code works than it is that test suite runs quickly. That said, it's nice when the test suite runs fast.

Plotting tests
~~~~~~~~~~~~~~

**TODO:** mainly this should cover `save_and_compare_images`
