from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).parent


def dogplot(*_, **__) -> None:
    """Shows who's a good boy"""
    pic = np.random.randint(1, 4)
    pic_path = HERE / "dogplot_images" / f"doggo_{pic}.jpg"

    img = plt.imread(pic_path)
    plt.imshow(img)
    plt.axis("off")
