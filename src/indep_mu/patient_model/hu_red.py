from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class CtCalibration:
    """Scanner-specific piecewise-linear HU -> RED calibration."""

    hu: np.ndarray
    red: np.ndarray

    def __post_init__(self) -> None:
        if self.hu.ndim != 1 or self.red.ndim != 1:
            raise ValueError("HU and RED arrays must be one-dimensional.")
        if self.hu.size != self.red.size:
            raise ValueError("HU and RED arrays must have equal length.")
        if self.hu.size < 2:
            raise ValueError("At least two calibration points are required.")
        if not np.all(np.diff(self.hu) > 0):
            raise ValueError("HU calibration points must be strictly increasing.")

    @property
    def hu_min(self) -> float:
        return float(self.hu[0])

    @property
    def hu_max(self) -> float:
        return float(self.hu[-1])

    def to_red(
        self,
        hu_values: np.ndarray,
        *,
        out_of_range: str = "raise",
    ) -> np.ndarray:
        """Convert HU to RED using linear interpolation.

        out_of_range='raise' rejects values outside the measured calibration
        interval. out_of_range='clip' reproduces Monaco pMC endpoint
        behaviour: HU below/above the CT-to-ED table receives the
        minimum/maximum RED respectively.

        The caller must choose clipping explicitly so the independent QA
        safety layer can still report out-of-range CT values.
        """
        values = np.asarray(hu_values, dtype=np.float64)

        if out_of_range not in {"raise", "clip"}:
            raise ValueError("out_of_range must be 'raise' or 'clip'.")

        below = values < self.hu_min
        above = values > self.hu_max

        if out_of_range == "raise" and (np.any(below) or np.any(above)):
            low = float(values.min())
            high = float(values.max())
            raise ValueError(
                f"HU outside calibration range [{self.hu_min}, {self.hu_max}]: "
                f"observed [{low}, {high}]"
            )

        if out_of_range == "clip":
            values = np.clip(values, self.hu_min, self.hu_max)

        return np.interp(values, self.hu, self.red)


DRT120KV = CtCalibration(
    hu=np.asarray(
        [-1000, -719, -66, -46, 3, 9, 11, 43, 45, 50, 55, 62, 119, 243, 277, 931, 1387, 2009],
        dtype=np.float64,
    ),
    red=np.asarray(
        [0.001, 0.227, 0.960, 0.979, 1.000, 1.007, 1.020, 1.041, 1.045, 1.047, 1.056, 1.057, 1.062, 1.119, 1.141, 1.469, 1.742, 2.335],
        dtype=np.float64,
    ),
)
