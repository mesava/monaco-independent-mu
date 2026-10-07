from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol, runtime_checkable

import numpy as np

from indep_mu.dicom.rtdose import RtDose

from .dose_quantity import DoseQuantity, require_matching_dose_quantity


@runtime_checkable
class PatientDoseSamplerGy(Protocol):
    frame_of_reference_uid: str
    dose_quantity: DoseQuantity

    def sample_patient_points_gy(self, points_patient_mm: np.ndarray) -> np.ndarray:
        ...


@dataclass(frozen=True)
class DoseDifferenceSummary:
    evaluated_voxel_count: int
    threshold_gy: float
    mean_difference_gy: float
    mean_absolute_difference_gy: float
    rms_difference_gy: float
    max_absolute_difference_gy: float
    mean_relative_difference_percent: float
    p95_absolute_relative_difference_percent: float


@dataclass(frozen=True)
class DoseDifferenceResult:
    dose_quantity: DoseQuantity
    reference_dose_gy: np.ndarray
    evaluated_dose_gy: np.ndarray
    difference_gy: np.ndarray
    relative_difference_percent: np.ndarray
    evaluation_mask: np.ndarray
    summary: DoseDifferenceSummary


def _dose_plane_patient_points(reference: RtDose, frame_index: int) -> np.ndarray:
    geometry = reference.geometry
    rows = geometry.rows
    columns = geometry.columns
    row_spacing, column_spacing = geometry.pixel_spacing_mm

    row_indices, column_indices = np.meshgrid(
        np.arange(rows, dtype=np.float64),
        np.arange(columns, dtype=np.float64),
        indexing="ij",
    )

    origin = geometry.frame_positions_patient_mm[frame_index]
    points = (
        origin[None, None, :]
        + column_indices[:, :, None]
        * column_spacing
        * geometry.column_direction[None, None, :]
        + row_indices[:, :, None]
        * row_spacing
        * geometry.row_direction[None, None, :]
    )
    return points.reshape((-1, 3))


def compare_on_rtdose_grid(
    reference: RtDose,
    evaluated: PatientDoseSamplerGy,
    *,
    reference_dose_quantity: DoseQuantity | None,
    threshold_fraction_of_reference_max: float = 0.10,
) -> DoseDifferenceResult:
    """Sample independent MC dose at native Monaco RTDOSE voxel centres."""

    dose_quantity = require_matching_dose_quantity(
        reference_dose_quantity,
        evaluated.dose_quantity,
    )

    if reference.frame_of_reference_uid != evaluated.frame_of_reference_uid:
        raise ValueError(
            "Reference RTDOSE and evaluated dose have different "
            "FrameOfReferenceUID values."
        )
    if not 0.0 <= threshold_fraction_of_reference_max < 1.0:
        raise ValueError(
            "threshold_fraction_of_reference_max must be in [0, 1)."
        )

    reference_gy = reference.dose_gy
    evaluated_gy = np.empty(reference_gy.shape, dtype=np.float64)

    for frame in range(reference.geometry.frames):
        points = _dose_plane_patient_points(reference, frame)
        values = np.asarray(
            evaluated.sample_patient_points_gy(points),
            dtype=np.float64,
        )
        expected_count = reference.geometry.rows * reference.geometry.columns
        if values.shape != (expected_count,):
            raise ValueError("Evaluated dose sampler returned unexpected shape.")
        evaluated_gy[frame] = values.reshape(
            (reference.geometry.rows, reference.geometry.columns)
        )

    difference = evaluated_gy - reference_gy
    threshold_gy = (
        float(np.max(reference_gy))
        * float(threshold_fraction_of_reference_max)
    )
    mask = reference_gy >= threshold_gy
    if not np.any(mask):
        raise ValueError("Dose threshold leaves no evaluation voxels.")

    relative = np.full(reference_gy.shape, np.nan, dtype=np.float64)
    safe = mask & (reference_gy > 0.0)
    relative[safe] = 100.0 * difference[safe] / reference_gy[safe]

    selected_difference = difference[mask]
    selected_relative = np.abs(relative[safe])

    summary = DoseDifferenceSummary(
        evaluated_voxel_count=int(np.count_nonzero(mask)),
        threshold_gy=threshold_gy,
        mean_difference_gy=float(np.mean(selected_difference)),
        mean_absolute_difference_gy=float(np.mean(np.abs(selected_difference))),
        rms_difference_gy=float(np.sqrt(np.mean(selected_difference**2))),
        max_absolute_difference_gy=float(np.max(np.abs(selected_difference))),
        mean_relative_difference_percent=float(np.mean(relative[safe])),
        p95_absolute_relative_difference_percent=float(
            np.percentile(selected_relative, 95)
        ),
    )

    return DoseDifferenceResult(
        dose_quantity=dose_quantity,
        reference_dose_gy=reference_gy.copy(),
        evaluated_dose_gy=evaluated_gy,
        difference_gy=difference,
        relative_difference_percent=relative,
        evaluation_mask=mask,
        summary=summary,
    )
