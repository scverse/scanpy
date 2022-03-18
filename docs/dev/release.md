# Making a release

This is a guide on how to make a release of scanpy for maintainers.

## Checking that PyPI will accept it

Sometimes PyPI will reject a build.
You can locally check whether the build is acceptable with `twine check`.

```shell
# First, make a build
python -m build

# Now check that build
twine check dist/*  # Assuming dist is otherwise empty
```

## Actually making release

First, make sure you're working repository is clean and is on the commit you'd like to release from.
Then follow these steps:

```shell
# Tag the commit with version info
git tag {version}  # where version is a version number like "1.7.0"

# Clear out old distributions
rm -r dist

# Build distributions and wheel
python -m build
# you can check that the previous step worked by installing from it's results
# e.g. pip install dist/scanpy-{version}-py3-none-any.whl

# Once you're confident the build looks good, push the tag to github
git push upstream {version}

# Upload wheel and code distribution to pypi
twine upload dist/scanpy-{version}
```

### Checking the distributions

If you're feeling cautious, you can:

- List the contents of a wheel with: `unzip -l dist/*.whl`
- Make a release candidate. Tag these with `{version}rc1` (e.g. `1.7.0rc1`)
- Upload them to <test.pypi.org> ([tutorial](https://packaging.python.org/en/latest/tutorials/packaging-projects/#uploading-the-distribution-archives))

## After making a release

After a major or minor release has been made:

- Tweet about it! Announce it on Zulip! Announce it on Discourse! Think about making a bot for this! Maybe actually do that?
- Create a new release notes file for the next minor release. This should only be added to the dev branch.
- Tag the development branch. If you just released `1.7.0`, this would be `1.8.0.dev0`.
- Create a new branch for this release series, like `1.7.x`. This should get a new release notes file.

After any release has been made:

- Create a new release notes file for the next bugfix release. This should be included in both dev and stable branches `release-latest.md`
- Create a milestone for the next release(s). For bugfix releases, this should have `on-merge: backport to 0.8.x` so [meeseeksdev](https://meeseeksbox.github.io) bot will create a backport PR. See {doc}`versioning` for more info.
