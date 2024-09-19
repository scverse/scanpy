(tests)=

# Tests

Possibly the most important part of contributing to any open source package is the test suite.
Implementations may change, but the only way we can know the code is working before making a release is the test suite.

## Running the tests

We use [pytest][] to test scanpy.
To run the tests, simply run `hatch test`.

It can take a while to run the whole test suite. There are a few ways to cut down on this while working on a PR:

1. Only run a subset of the tests.
   This can be done by specifying paths or test name patterns using the `-k` argument (e.g. `hatch test test_plotting.py` or `hatch test -k "test_umap*"`)
2. Run the tests in parallel using the `-n` argument (e.g. `hatch test -n 8`).

[pytest]: https://docs.pytest.org/en/stable/

### Miscellaneous tips

- A lot of warnings can be thrown while running the test suite.
  It's often easier to read the test results with them hidden via the `--disable-pytest-warnings` argument.

## Writing tests

You can refer to the [existing test suite][] for examples.
If you haven't written tests before, Software Carpentry has an [in-depth testing guide][].

We highly recommend using [Test-Driven Development][] when contributing code.
This not only ensures you have tests written, it often makes implementation easier since you start out with a specification for your function.

Consider parameterizing your tests using the `pytest.mark.parameterize` and `pytest.fixture` decorators.
You can read more about [fixtures][] in pytest’s documentation, but we’d also recommend searching our test suite for existing usage.

[existing test suite]: https://github.com/scverse/scanpy/tree/main/scanpy/tests
[in-depth testing guide]: https://katyhuff.github.io/2016-07-11-scipy/testing/
[test-driven development]: https://en.wikipedia.org/wiki/Test-driven_development
[fixtures]: https://docs.pytest.org/en/stable/fixture.html

### What to test

If you're not sure what to tests about your function, some ideas include:

- Are there arguments which conflict with each other? Check that if they are both passed, the function throws an error (see [`pytest.raises`][] docs).
- Are there input values which should cause your function to error?
- Did you add a helpful error message that recommends better outputs? Check that that error message is actually thrown.
- Can you place bounds on the values returned by your function?
- Are there different input values which should generate equivalent output (e.g. if an array is sparse or dense)?
- Do you have arguments which should have orthogonal effects on the output? Check that they are independent. For example, if there is a flag for extended output, the base output should remain the same either way.
- Are you optimizing a method? Check that it's results are the same as a gold standard implementation.

[`pytest.raises`]: https://docs.pytest.org/en/stable/assert.html#assertions-about-expected-exceptions

### Performance

It's more important that you're accurately testing the code works than it is that test suite runs quickly.
That said, it's nice when the test suite runs fast.

You can check how long tests take to run by passing `--durations=0` argument to `pytest`.
Hopefully your new tests won't show up on top!
Some approaches to this include:

- Is there a common setup/ computation happening in each test? Consider caching these in a [scoped test fixture][].
- Is the behaviour you're testing for dependent on the size of the data? If not, consider reducing it.

[scoped test fixture]: https://docs.pytest.org/en/stable/fixture.html#sharing-test-data

### Plotting tests

While computational functions will return arrays and values, it can be harder to work with the output of plotting functions.

To make this easier, we use the `image_comparer` fixture for comparing plotting results (search the test suite for example usage).
This is used to check that generated plots look the same as they did previously.
Reference images (the expected output) are stored as `expected.png` to relevant tests directory under `scanpy/tests/_images`.
When run, the test suite will generate `actual.png` files for each check.
These files are compared, and if the `actual` plot differs from the reference plot, a `diff` of the images is also generated.
Paths for all these files will be reported when a test fails, and images for failed plots can be viewed via the :doc:`CI interface <ci>`.

A common gotcha here is that plots often change slightly on different machines/ OSs.
`scanpy`'s test suite sets a number of environment variables to ensure as similar of plots as possible.
When adding new reference plots, the recommended workflow is to write the test as though an expected result already exists, run it once to generate the output, then move that output to the reference directory.
