from __future__ import annotations

from importlib.resources import files

import matplotlib.pyplot as plt
import numpy as np


def dogplot(*_, **__) -> None:
    """Show who's a good boy."""
    pic = np.random.randint(1, 4)
    img_path = files("scanpy.plotting") / f"dogplot_images/doggo_{pic}.webp"

    plt.figure(figsize=(3, 3))

    with img_path.open("rb") as f:
        img = plt.imread(f)
    plt.imshow(img)

    plt.axis("off")
