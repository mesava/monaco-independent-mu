from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class DicomBankPositions:
    """Two DICOM leaf/jaw banks in IEC subscript order.

    DICOM stores 2N values as 101..1N followed by 201..2N.  This class keeps
    that bank identity intact; it does not guess which bank is the negative
    side of an opening.
    """

    bank1_mm: np.ndarray
    bank2_mm: np.ndarray

    def __post_init__(self) -> None:
        bank1 = np.asarray(self.bank1_mm, dtype=np.float64)
        bank2 = np.asarray(self.bank2_mm, dtype=np.float64)
        if bank1.ndim != 1 or bank2.ndim != 1:
            raise ValueError("DICOM bank arrays must be one-dimensional.")
        if bank1.shape != bank2.shape:
            raise ValueError("DICOM bank arrays must have equal length.")
        if bank1.size == 0:
            raise ValueError("DICOM bank arrays must not be empty.")
        if not np.all(np.isfinite(bank1)) or not np.all(np.isfinite(bank2)):
            raise ValueError("DICOM bank positions must be finite.")


def split_dicom_banks(
    positions_mm: tuple[float, ...] | list[float] | np.ndarray,
    *,
    pair_count: int,
) -> DicomBankPositions:
    """Split LeafJawPositions into DICOM bank 1 and bank 2.

    The DICOM order is 101,102,...,1N,201,202,...,2N.
    """

    values = np.asarray(positions_mm, dtype=np.float64)
    if values.ndim != 1:
        raise ValueError("LeafJawPositions must be one-dimensional.")
    if pair_count < 1:
        raise ValueError("pair_count must be positive.")
    if values.size != 2 * pair_count:
        raise ValueError(
            f"Expected {2 * pair_count} LeafJawPositions, got {values.size}."
        )

    return DicomBankPositions(
        bank1_mm=values[:pair_count].copy(),
        bank2_mm=values[pair_count:].copy(),
    )


@dataclass(frozen=True)
class IsocenterProjection:
    """Project an IEC opening coordinate from isocenter to a source plane.

    This helper is deliberately explicit: callers must first establish that the
    DICOM opening coordinate they supply is an isocenter-projected coordinate.

    The projection assumes the radiation source/focus is at z=0 and the
    isocenter is at SAD.
    """

    sad_mm: float

    def __post_init__(self) -> None:
        if self.sad_mm <= 0:
            raise ValueError("SAD must be positive.")

    def to_source_distance_cm(
        self,
        position_mm: np.ndarray | float,
        *,
        source_distance_cm: float,
    ) -> np.ndarray:
        if source_distance_cm <= 0:
            raise ValueError("source_distance_cm must be positive.")
        position = np.asarray(position_mm, dtype=np.float64)
        if not np.all(np.isfinite(position)):
            raise ValueError("Position values must be finite.")
        scale = (source_distance_cm * 10.0) / self.sad_mm
        return (position * scale) / 10.0

    def to_isocenter_mm(
        self,
        position_cm: np.ndarray | float,
        *,
        source_distance_cm: float,
    ) -> np.ndarray:
        if source_distance_cm <= 0:
            raise ValueError("source_distance_cm must be positive.")
        position = np.asarray(position_cm, dtype=np.float64)
        if not np.all(np.isfinite(position)):
            raise ValueError("Position values must be finite.")
        scale = self.sad_mm / (source_distance_cm * 10.0)
        return position * 10.0 * scale
