#!/usr/bin/env python3
"""Write a symmetric Bethe hybridization on the reduced broadening mesh."""

import argparse
import math
from pathlib import Path


def mesh_points():
    magnitudes = [1e-4]
    while magnitudes[-1] * 1.05 < 40.0:
        magnitudes.append(magnitudes[-1] * 1.05)
    magnitudes.extend((1.0, 40.0))
    positive = sorted(set(magnitudes))
    return [-value for value in reversed(positive)] + positive


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("gamma", type=Path)
    parser.add_argument("mesh", type=Path)
    arguments = parser.parse_args()

    points = mesh_points()
    gamma_rows = []
    mesh_rows = []
    for omega in points:
        gamma = (
            0.5 * math.sqrt(max(0.0, 1.0 - omega * omega))
            if abs(omega) < 1.0
            else 0.0
        )
        gamma_rows.append(f"{omega:.17g} {gamma:.17g}\n")
        mesh_rows.append(f"{omega:.17g} 0\n")
    arguments.gamma.write_text("".join(gamma_rows), encoding="ascii")
    arguments.mesh.write_text("".join(mesh_rows), encoding="ascii")


if __name__ == "__main__":
    main()
