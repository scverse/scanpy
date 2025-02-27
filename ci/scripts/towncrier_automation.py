#!/usr/bin/env python3
# /// script
# dependencies = [ "towncrier", "packaging" ]
# ///
"""Script to automate towncrier release note PRs."""

from __future__ import annotations

import argparse
import subprocess
from typing import TYPE_CHECKING

from packaging.version import Version

if TYPE_CHECKING:
    from collections.abc import Sequence


class Args(argparse.Namespace):
    """Command line arguments."""

    version: str
    dry_run: bool


def parse_args(argv: Sequence[str] | None = None) -> Args:
    """Construct a CLI argument parser."""
    parser = argparse.ArgumentParser(
        prog="towncrier-automation",
        description=(
            "This script runs towncrier for a given version, "
            "creates a branch off of the current one, "
            "and then creates a PR into the original branch with the changes. "
            "The PR will be backported to main if the current branch is not main."
        ),
    )
    parser.add_argument(
        "version",
        type=str,
        help=(
            "The new version for the release must have at least three parts, like `major.minor.patch` and no `major.minor`. "
            "It can have a suffix like `major.minor.patch.dev0` or `major.minor.0rc1`."
        ),
    )
    parser.add_argument(
        "--dry-run",
        help="Whether or not to dry-run the actual creation of the pull request",
        action="store_true",
    )
    args = parser.parse_args(argv, Args())
    # validate the version
    if len(Version(args.version).release) != 3:
        msg = f"Version argument {args.version} must contain major, minor, and patch version."
        raise ValueError(msg)
    return args


def main(argv: Sequence[str] | None = None) -> None:
    """Run main entry point."""
    args = parse_args(argv)

    # Run towncrier
    subprocess.run(
        ["towncrier", "build", f"--version={args.version}", "--yes"], check=True
    )

    # Check if we are on the main branch to know if we need to backport
    base_branch = subprocess.run(
        ["git", "rev-parse", "--abbrev-ref", "HEAD"],
        capture_output=True,
        text=True,
        check=True,
    ).stdout.strip()
    pr_description = "" if base_branch == "main" else "@meeseeksdev backport to main"
    branch_name = f"release_notes_{args.version}"

    # Create a new branch + commit
    subprocess.run(["git", "switch", "-c", branch_name], check=True)
    subprocess.run(["git", "add", "docs/release-notes"], check=True)
    pr_title = f"(chore): generate {args.version} release notes"
    subprocess.run(["git", "commit", "-m", pr_title], check=True)

    # push
    if not args.dry_run:
        subprocess.run(
            ["git", "push", "--set-upstream", "origin", branch_name], check=True
        )
    else:
        print("Dry run, not pushing")

    # Create a PR
    subprocess.run(
        [
            "gh",
            "pr",
            "create",
            f"--base={base_branch}",
            f"--title={pr_title}",
            f"--body={pr_description}",
            *(
                ["--label=no milestone", "--label=Development Process 🚀"]
                if base_branch == "main"
                else []
            ),
            *(["--dry-run"] if args.dry_run else []),
        ],
        check=True,
    )

    # Enable auto-merge
    if not args.dry_run:
        subprocess.run(
            ["gh", "pr", "merge", branch_name, "--auto", "--squash"], check=True
        )
    else:
        print("Dry run, not merging")


if __name__ == "__main__":
    main()
