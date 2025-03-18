"""Test our test helpers."""

from __future__ import annotations

import numpy as np

from testing.scanpy._helpers import random_mask


def test_random_mask():
    ns_true = np.array([int(random_mask(4).sum()) for _ in range(1000)])
    np.testing.assert_equal(ns_true, [2] * 1000)
