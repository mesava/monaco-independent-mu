from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class RectilinearDoseGrid:
    """Dose on a rectilinear grid of voxel-centre coordinates.

    Internal array convention is dose[z, y, x].
    Coordinate units are cm; dose units are supplied by the caller and are
    preserved by interpolation.
    """

    x_cm: np.ndarray
    y_cm: np.ndarray
    z_cm: np.ndarray
    dose: np.ndarray

    def __post_init__(self) -> None:
        x = np.asarray(self.x_cm, dtype=np.float64)
        y = np.asarray(self.y_cm, dtype=np.float64)
        z = np.asarray(self.z_cm, dtype=np.float64)
        dose = np.asarray(self.dose, dtype=np.float64)

        for name, vector in (("x", x), ("y", y), ("z", z)):
            if vector.ndim != 1 or vector.size < 2:
                raise ValueError(f"{name} coordinates must be a 1-D array with >=2 points.")
            if np.any(np.diff(vector) <= 0):
                raise ValueError(f"{name} coordinates must be strictly increasing.")

        expected = (z.size, y.size, x.size)
        if dose.shape != expected:
            raise ValueError(
                f"Dose shape {dose.shape} does not match grid shape {expected}."
            )
        if not np.all(np.isfinite(dose)):
            raise ValueError("Dose grid contains non-finite values.")


def _bracket(
    coordinates: np.ndarray,
    values: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return lower/upper indices and interpolation fraction."""

    if np.any(values < coordinates[0]) or np.any(values > coordinates[-1]):
        raise ValueError("Interpolation point lies outside dose-grid centres.")

    upper = np.searchsorted(coordinates, values, side="right")
    upper = np.clip(upper, 1, coordinates.size - 1)
    lower = upper - 1

    denominator = coordinates[upper] - coordinates[lower]
    fraction = (values - coordinates[lower]) / denominator
    return lower, upper, fraction


def trilinear_sample(
    grid: RectilinearDoseGrid,
    points_cm: np.ndarray,
) -> np.ndarray:
    """Sample dose at N points using trilinear interpolation."""

    points = np.asarray(points_cm, dtype=np.float64)
    if points.ndim != 2 or points.shape[1] != 3:
        raise ValueError("points_cm must have shape (N, 3) ordered as x,y,z.")
    if not np.all(np.isfinite(points)):
        raise ValueError("Interpolation points must be finite.")

    x0, x1, fx = _bracket(np.asarray(grid.x_cm), points[:, 0])
    y0, y1, fy = _bracket(np.asarray(grid.y_cm), points[:, 1])
    z0, z1, fz = _bracket(np.asarray(grid.z_cm), points[:, 2])

    dose = np.asarray(grid.dose, dtype=np.float64)

    c000 = dose[z0, y0, x0]
    c100 = dose[z0, y0, x1]
    c010 = dose[z0, y1, x0]
    c110 = dose[z0, y1, x1]
    c001 = dose[z1, y0, x0]
    c101 = dose[z1, y0, x1]
    c011 = dose[z1, y1, x0]
    c111 = dose[z1, y1, x1]

    c00 = c000 * (1.0 - fx) + c100 * fx
    c10 = c010 * (1.0 - fx) + c110 * fx
    c01 = c001 * (1.0 - fx) + c101 * fx
    c11 = c011 * (1.0 - fx) + c111 * fx

    c0 = c00 * (1.0 - fy) + c10 * fy
    c1 = c01 * (1.0 - fy) + c11 * fy
    return c0 * (1.0 - fz) + c1 * fz
