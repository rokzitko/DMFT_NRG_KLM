"""Shared edge-resolved Bethe-lattice DOS and transport tables."""

import math
import os


PHYSICAL_POINTS = 2001
PADDING_START = 1300
PADDING_EDGE = 1000


def publication_mode(path):
    """Preserve an output's mode, or apply the process umask for a new file."""
    try:
        return path.stat().st_mode & 0o777
    except FileNotFoundError:
        mask = os.umask(0)
        os.umask(mask)
        return 0o666 & ~mask


def bethe_rows():
    """Return matching (epsilon, rho, phi) rows on a theta-refined mesh."""
    rows = [
        (index / 1000, 0.0, 0.0)
        for index in range(-PADDING_START, -PADDING_EDGE)
    ]

    denominator = PHYSICAL_POINTS - 1
    for index in range(PHYSICAL_POINTS):
        theta = -math.pi / 2 + math.pi * index / denominator
        epsilon = math.sin(theta)
        cosine = max(0.0, math.cos(theta))
        if index == 0:
            epsilon = -1.0
            cosine = 0.0
        elif index == denominator:
            epsilon = 1.0
            cosine = 0.0
        rho = 2.0 / math.pi * cosine
        rows.append((epsilon, rho, rho * cosine * cosine))

    rows.extend(
        (index / 1000, 0.0, 0.0)
        for index in range(PADDING_EDGE + 1, PADDING_START + 1)
    )

    if len(rows) != 2601:
        raise RuntimeError(f"Unexpected Bethe table size: {len(rows)}")
    if any(first[0] >= second[0] for first, second in zip(rows, rows[1:])):
        raise RuntimeError("Bethe table mesh is not strictly increasing")
    return rows
