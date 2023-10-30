# Making a release

1. Go to GitHub’s [releases][] page
2. Click “Draft a new release”
3. Click “Choose a tag” and type the version of the tag you want to release, such as “1.9.6”
4. Click “**+ Create new tag: 1.x.x** on publish”
5. If the version is a pre-release version, such as “1.10.0a1”, tick the “Set as a pre-release” checkbox

[releases]: https://github.com/scverse/scanpy/releases

## After making a release

After a major or minor release has been made:

- Tweet about it! Announce it on Zulip! Announce it on Discourse! Think about making a bot for this! Maybe actually do that?
- Create a new release notes file for the next minor release. This should only be added to the dev branch.
- Tag the development branch. If you just released `1.7.0`, this would be `1.8.0.dev0`.
- Create a new branch for this release series, like `1.7.x`. This should get a new release notes file.

After any release has been made:

- Create a new release notes file for the next bugfix release.
  This should be included in both dev and stable branches.
- Create a milestone for the next release(s).
  For bugfix releases, this should have `on-merge: backport to 0.8.x`,
  so the [meeseeksdev][] bot will create a backport PR. See {doc}`versioning` for more info.
- Clear out and close the milestone you just made a release for.

[meeseeksdev]: https://meeseeksbox.github.io
