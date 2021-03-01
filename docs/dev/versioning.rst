Versioning
==========

.. note::

    We are currently experimenting with our development practices.
    These are currently documented on a best effort basis, but may not be completely accurate.

We try to follow `semantic versioning <https://semver.org>`__ with our versioning scheme.
This scheme breaks down a version number into `{major.minor.point}` sections.
At a `point` release, there should be no changes beyond bug fixes.
`minor` releases can include new features.
`major` releases can break old APIs.

To be sure we can follow this scheme and maintain some agility in development, we use some tooling and development practices.
When a minor release is made, a release branch should be cut and pushed to the main repo (e.g. `1.7.x` for the `1.7` release series).

For PRs which fix an bug in the most recent minor release, the changes will need to added to both the development and release branches.
To accomplish this, PRs which fix bugs must be labelled as such.
After approval, a developer will notify the `meeseeks bot <https://meeseeksbox.github.io>`__ to open a backport PR onto the release branch via a comment saying:

    @Meeseeksdev backport <branch>

Where "<branch>" is the most recent release branch.

The bot will attempt to make a backport and open a PR.
This will sometimes require manual intervention due to merge conflicts or test failures.
