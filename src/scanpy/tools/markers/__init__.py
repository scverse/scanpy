"""Rank marker genes that distinguish groups of cells."""

from __future__ import annotations

from ._api import logreg, ttest, wilcoxon

__all__ = ["logreg", "ttest", "wilcoxon"]
