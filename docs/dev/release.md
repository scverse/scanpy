# Making a release

First, check out {doc}`versioning` to see which kind of release you want to make.
That page also explains concepts like *pre-releases* and applications thereof.

## Actually making the release

1. Go to GitHub’s [releases][] page
2. Click “Draft a new release”
3. Click “Choose a tag” and type the version of the tag you want to release, such as `1.9.6`
4. Click “**+ Create new tag: 1.\<minor>.\<patch>** on publish”
5. If the version is a *pre-release* version, such as `1.7.0rc1` or `1.10.0a1`, tick the “Set as a pre-release” checkbox

[releases]: https://github.com/scverse/scanpy/releases

## After making a release

After *any* release has been made:

- Create a new release notes file for the next bugfix release.
  This should be included in both dev and stable branches.
- Create a milestone for the next release (in case you made a bugfix release) or releases (in case of a major/minor release).
  For bugfix releases, this should have `on-merge: backport to 0.<minor>.x`,
  so the [meeseeksdev][] bot will create a backport PR. See {doc}`versioning` for more info.
- Clear out and close the milestone you just made a release for.

After a *major* or *minor* release has been made:

- Tweet about it! Announce it on Zulip! Announce it on Discourse! Think about making a bot for this! Maybe actually do that?
- Create a new release notes file for the next minor release. This should only be added to the dev branch.
- Tag the development branch. If you just released `1.7.0`, this would be `1.8.0.dev0`.
- Create a new branch for this release series, like `1.7.x`. This should get a new release notes file.

[meeseeksdev]: https://meeseeksbox.github.io

## Debugging the build process

If you changed something about the build process (e.g. [Hatchling’s build configuration][hatch-build]),
or something about the package’s structure,
you might want to manually check if the build and upload process behaves as expected:

```shell
# Clear out old distributions
rm -r dist

# Build source distribution and wheel both
python -m build

# Now check those build artifacts
twine check dist/*

# List the wheel archive’s contents
bsdtar -tf dist/*.whl

```

You can also upload the package to <test.pypi.org> ([tutorial][testpypi tutorial])

[testpypi tutorial]: https://packaging.python.org/en/latest/tutorials/packaging-projects/#uploading-the-distribution-archives

```
twine upload --repository testpypi dist/*
```

The above approximates what the [publish workflow][] does automatically for us.
If you want to replicate the process more exactly, make sure you are careful,
and create a version tag before building (make sure you delete it after uploading to TestPyPI!).

[hatch-build]: https://hatch.pypa.io/latest/config/build/
[publish workflow]: https://github.com/scverse/scanpy/tree/main/.github/workflows/publish.yml
