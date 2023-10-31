# Versioning

```{note}
We are currently experimenting with our development practices.
These are currently documented on a best effort basis, but may not be completely accurate.
```

## Semantic versioning

We try to follow [semantic versioning](https://semver.org) with our versioning scheme.
This scheme breaks down a version number into `{major.minor.point}` sections.
At a `point` release, there should be no changes beyond bug fixes.
`minor` releases can include new features.
`major` releases can break old APIs.

### Version numbers

Valid version numbers are described in [PEP 440](https://peps.python.org/pep-0440/).

[Pre-releases](https://peps.python.org/pep-0440/#pre-releases)
:   should have versions like `1.7.0rc1` or `1.7.0rc2`.
[Development versions](https://peps.python.org/pep-0440/#developmental-releases)
:   should look like `1.8.0.dev0`, with a commit hash optionally appended as a local version identifier (e.g. `1.8.0.dev2+g00ad77b`).

## Tooling

To be sure we can follow this scheme and maintain some agility in development, we use some tooling and development practices.
When a minor release is made, a release branch should be cut and pushed to the main repo (e.g. `1.7.x` for the `1.7` release series).

For PRs which fix an bug in the most recent minor release, the changes will need to added to both the development and release branches.
To accomplish this, PRs which fix bugs must be labelled as such.
After approval, a developer will notify the [meeseeks bot](https://meeseeksbox.github.io) to open a backport PR onto the release branch via a comment saying:

> @Meeseeksdev backport \<branch>

Where "\<branch>" is the most recent release branch.

The bot will attempt to make a backport and open a PR.
This will sometimes require manual intervention due to merge conflicts or test failures.
