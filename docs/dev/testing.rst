.. _tests:

Tests
=====

Possibly the most important part of contributing to any open source package is the test suite.
Implementations may change, but the only way we can know the code is working before making a release is the test suite.

Running the tests
-----------------

We use `pytest <https://docs.pytest.org/en/stable/>`__ to test scanpy.
To run the tests first make sure you have the required dependencies (``pip install -e ".[tests]"``), then run ``pytest`` from the root of the repository.

It can take a while to run the whole test suite. There are a few ways to cut down on this while working on a PR:

1. Only run a subset of the tests. This can be done with the `-k` argument from pytest (e.g. ``pytest -k test_plotting.py`` or ``pytest -k "test_umap*"``
2. Run the tests in parallel. If you install the pytest extension `pytest-xdist <https://github.com/pytest-dev/pytest-xdist>`__ you can run tests in parallel with the ``--numprocesses`` argument to pytest (e.g. ``pytest -n 8``).

Miscellaneous tips
~~~~~~~~~~~~~~~~~~

- A lot of warnings can be thrown while running the test suite. It's often easier to read the test results with them hidden via the `--disable-pytest-warnings` argument.

Writing tests
-------------

You can refer to the `existing test suite <https://github.com/theislab/scanpy/tree/master/scanpy/tests>`__ for examples.
If you haven't written tests before, Software Carpentry has an `in-depth guide <http://katyhuff.github.io/python-testing/>`__ on the topic.

We highly recommend using `Test Driven Development <https://en.wikipedia.org/wiki/Test-driven_development>`__ when contributing code.
This not only ensures you have tests written, it often makes implementation easier since you start out with a specification for your function.

Consider parameterizing your tests using the `pytest.mark.parameterize` and `pytest.fixture` decorators.
Documentation on these can be found `here <https://docs.pytest.org/en/stable/fixture.html>`__, but we'd also recommend searching our test suite for existing usage.

What to test
~~~~~~~~~~~~

If you're not sure what to tests about your function, some ideas include:

- Are there arguments which conflict with each other? Check that if they are both passed, the function throws an error (see `pytest.raises` `in the pytest docs <https://docs.pytest.org/en/stable/assert.html#assertions-about-expected-exceptions>`__.
- Are there input values which should cause your function to error?
- Did you add a helpful error message that recommends better outputs? Check that that error message is actually thrown.
- Can you place bounds on the values returned by your function?
- Are there different input values which should generate equivalent output (e.g. if an array is sparse or dense)?
- Do you have arguments which should have orthogonal effects on the output? Check that they are independent. For example, if there is a flag for extended output, the base output should remain the same either way.
- Are you optimizing a method? Check that it's results are the same as a gold standard implementation.

Performance
~~~~~~~~~~~

It's more important that you're accurately testing the code works than it is that test suite runs quickly.
That said, it's nice when the test suite runs fast.

You can check how long tests take to run by passing `--durations=0` argument to `pytest`.
Hopefully your new tests won't show up on top!
Some approaches to this include:

- Is there a common setup/ computation happening in each test? Consider caching these in a `scoped test fixture <https://docs.pytest.org/en/stable/fixture.html#sharing-test-data>`__.
- Is the behaviour you're testing for dependent on the size of the data? If not, consider reducing it.

Plotting tests
~~~~~~~~~~~~~~

While computational functions will return arrays and values, it can be harder to work with the output of plotting functions.

To make this easier, we use the `image_comparer` fixture for comparing plotting results (search the test suite for example usage).
This is used to check that generated plots look the same as they did previously.
Reference images (the expected output) are kept `scanpy/tests/_images` and compared with the results of running the test suite by running the tests.

A common gotcha here is that plots often change slightly on different machines/ OSs.
`scanpy`'s test suite sets a number of environment variables to ensure as similar of plots as possible.
When adding new reference plots, the recommended workflow is to write the test as though an expected result already exists, run it once to generate the output, then move that output to the reference directory.
