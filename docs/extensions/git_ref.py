from __future__ import annotations

from functools import lru_cache

import re
import subprocess
from sphinx.application import Sphinx


def git(*args: str) -> str:
    return subprocess.check_output(["git", *args]).strip().decode()


# https://github.com/DisnakeDev/disnake/blob/7853da70b13fcd2978c39c0b7efa59b34d298186/docs/conf.py#L192
@lru_cache
def get() -> str | None:
    """Current git reference. Uses branch/tag name if found, otherwise uses commit hash"""
    git_ref = None
    try:
        git_ref = git("name-rev", "--name-only", "--no-undefined", "HEAD")
        git_ref = re.sub(r"^(remotes/[^/]+|tags)/", "", git_ref)
    except Exception:
        pass

    # (if no name found or relative ref, use commit hash instead)
    if not git_ref or re.search(r"[\^~]", git_ref):
        try:
            git_ref = git("rev-parse", "HEAD")
        except Exception:
            git_ref = "master"
    return git_ref


def set_ref(app: Sphinx, pagename, templatename, context, doctree):
    context['theme_repository_branch'] = get()


def setup(app: Sphinx) -> None:
    app.connect("html-page-context", set_ref)
