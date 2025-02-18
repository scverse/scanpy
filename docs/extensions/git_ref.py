"""Extension to inject ``html_theme_options["repository_branch"]``."""

from __future__ import annotations

import re
import subprocess
from functools import lru_cache
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from sphinx.application import Sphinx
    from sphinx.config import Config


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
            git_ref = "main"
    return git_ref


def set_ref(app: Sphinx, config: Config):
    app.config["html_theme_options"]["repository_branch"] = get() or "main"


def setup(app: Sphinx) -> None:
    app.connect("config-inited", set_ref)
