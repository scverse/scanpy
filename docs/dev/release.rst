Making a release
================

This is a guide on how to make a release of scanpy for maintainers.

First, make sure you're working repository is clean and is on the commit you'd like to release from.
Then follow these steps:

.. code:: shell

    # Tag the commit with version info
    git tag {version}  # where version is a version number like "1.7.0"

    # Build distributions and wheel
    python setup.py sdist bdist_wheel
    # you can check that the previous step worked by installing from it's results
    # e.g. pip install dist/scanpy-{version}-py3-none-any.whl

    # Once you're confident the build looks good, push the tag to github
    git push upstream {version}

    # Upload wheel and code distribution to pypi
    twine upload dist/scanpy-{version}
