from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass(frozen=True)
class ThreeDDose:
    """Dense ASCII DOSXYZnrc .3ddose data.

    Arrays use project convention (z, y, x); x is the fastest-varying value in
    the file, matching DOSXYZnrc output order.

    The numeric dose is kept in the native DOSXYZnrc normalization.  This class
    deliberately does not label it as an absolute clinical Gy/MU value until an
    independent calibration factor has been applied.
    """

    x_bound_cm: np.ndarray
    y_bound_cm: np.ndarray
    z_bound_cm: np.ndarray
    dose: np.ndarray
    relative_uncertainty: np.ndarray

    def __post_init__(self) -> None:
        x = np.asarray(self.x_bound_cm, dtype=np.float64)
        y = np.asarray(self.y_bound_cm, dtype=np.float64)
        z = np.asarray(self.z_bound_cm, dtype=np.float64)
        dose = np.asarray(self.dose, dtype=np.float64)
        uncertainty = np.asarray(self.relative_uncertainty, dtype=np.float64)

        if any(vector.ndim != 1 for vector in (x, y, z)):
            raise ValueError("3ddose boundary arrays must be one-dimensional.")
        if any(vector.size < 2 for vector in (x, y, z)):
            raise ValueError("Each 3ddose axis requires at least two boundaries.")
        if any(np.any(np.diff(vector) <= 0) for vector in (x, y, z)):
            raise ValueError("3ddose boundaries must be strictly increasing.")

        expected = (z.size - 1, y.size - 1, x.size - 1)
        if dose.shape != expected:
            raise ValueError(
                f"Dose shape {dose.shape} does not match boundaries {expected}."
            )
        if uncertainty.shape != expected:
            raise ValueError(
                "Relative-uncertainty shape does not match dose geometry."
            )
        if not np.all(np.isfinite(dose)):
            raise ValueError("3ddose dose array contains non-finite values.")
        if not np.all(np.isfinite(uncertainty)):
            raise ValueError("3ddose uncertainty array contains non-finite values.")
        if np.any(uncertainty < 0):
            raise ValueError("3ddose relative uncertainties must be non-negative.")

    @property
    def shape(self) -> tuple[int, int, int]:
        return tuple(int(value) for value in self.dose.shape)

    @property
    def voxel_centres_cm(
        self,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        return (
            0.5 * (self.x_bound_cm[:-1] + self.x_bound_cm[1:]),
            0.5 * (self.y_bound_cm[:-1] + self.y_bound_cm[1:]),
            0.5 * (self.z_bound_cm[:-1] + self.z_bound_cm[1:]),
        )

    def scaled(self, factor: float) -> "ThreeDDose":
        """Scale dose only; relative statistical uncertainty is unchanged."""

        if not np.isfinite(factor) or factor <= 0:
            raise ValueError("Dose scaling factor must be finite and positive.")
        return ThreeDDose(
            x_bound_cm=np.asarray(self.x_bound_cm, dtype=np.float64).copy(),
            y_bound_cm=np.asarray(self.y_bound_cm, dtype=np.float64).copy(),
            z_bound_cm=np.asarray(self.z_bound_cm, dtype=np.float64).copy(),
            dose=np.asarray(self.dose, dtype=np.float64) * factor,
            relative_uncertainty=np.asarray(
                self.relative_uncertainty,
                dtype=np.float64,
            ).copy(),
        )


def _tokens(path: str | Path) -> list[str]:
    return Path(path).read_text(encoding="ascii").split()


def read_3ddose(path: str | Path) -> ThreeDDose:
    """Read the standard dense ASCII DOSXYZnrc .3ddose format."""

    tokens = _tokens(path)
    if len(tokens) < 3:
        raise ValueError("3ddose file is incomplete.")

    position = 0

    def take_int() -> int:
        nonlocal position
        if position >= len(tokens):
            raise ValueError("Unexpected end of 3ddose file.")
        value = int(tokens[position])
        position += 1
        return value

    def take_float(count: int) -> np.ndarray:
        nonlocal position
        stop = position + count
        if stop > len(tokens):
            raise ValueError("Unexpected end of 3ddose file.")
        values = np.asarray(
            [float(value) for value in tokens[position:stop]],
            dtype=np.float64,
        )
        position = stop
        return values

    nx, ny, nz = take_int(), take_int(), take_int()
    if nx < 1 or ny < 1 or nz < 1:
        raise ValueError("3ddose voxel dimensions must be positive.")

    x_bound = take_float(nx + 1)
    y_bound = take_float(ny + 1)
    z_bound = take_float(nz + 1)

    count = nx * ny * nz
    dose = take_float(count).reshape((nz, ny, nx))
    uncertainty = take_float(count).reshape((nz, ny, nx))

    if position != len(tokens):
        raise ValueError(
            f"Unexpected trailing values in 3ddose file: {len(tokens) - position}."
        )

    return ThreeDDose(
        x_bound_cm=x_bound,
        y_bound_cm=y_bound,
        z_bound_cm=z_bound,
        dose=dose,
        relative_uncertainty=uncertainty,
    )


def write_3ddose(path: str | Path, data: ThreeDDose) -> None:
    """Write a deterministic standard dense ASCII .3ddose file."""

    nz, ny, nx = data.shape
    lines = [f"{nx} {ny} {nz}"]

    for vector in (data.x_bound_cm, data.y_bound_cm, data.z_bound_cm):
        lines.append(" ".join(f"{float(value):.12g}" for value in vector))

    lines.append(
        " ".join(
            f"{float(value):.12g}"
            for value in np.asarray(data.dose).reshape(-1)
        )
    )
    lines.append(
        " ".join(
            f"{float(value):.12g}"
            for value in np.asarray(data.relative_uncertainty).reshape(-1)
        )
    )

    Path(path).write_text("\n".join(lines) + "\n", encoding="ascii")
