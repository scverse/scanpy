from __future__ import annotations

from importlib.resources import files

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.image import imread


def dogplot(n: int | None = None, *, show: bool = False) -> None:
    """Show who’s a good boy."""
    rng = np.random.default_rng()
    n = int(rng.integers(1, 4))
    img_path = files("scanpy.plotting") / f"dogplot_images/doggo_{n}.webp"
    with img_path.open("rb") as f:
        img = imread(f)

    _, ax = plt.subplots(figsize=(3, 3))
    ax.imshow(img)
    ax.set_axis_off()
