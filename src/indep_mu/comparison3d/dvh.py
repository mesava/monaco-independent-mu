from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol, runtime_checkable

import numpy as np

from indep_mu.dicom.ct import CtSeries


@runtime_checkable
class PatientDoseSamplerGy(Protocol):
    frame_of_reference_uid: str

    def sample_patient_points_gy(self, points_patient_mm: np.ndarray) -> np.ndarray:
        ...


@dataclass(frozen=True)
class StructureDoseSamples:
    dose_gy: np.ndarray
    voxel_volume_cm3: float

    @property
    def volume_cm3(self) -> float:
        return float(self.dose_gy.size * self.voxel_volume_cm3)


@dataclass(frozen=True)
class DvhMetrics:
    volume_cm3: float
    d98_gy: float
    d95_gy: float
    d50_gy: float
    d2_gy: float
    mean_gy: float
    min_gy: float
    max_gy: float


def _ct_mask_patient_points(ct: CtSeries, mask: np.ndarray) -> np.ndarray:
    selected = np.asarray(mask, dtype=bool)
    if selected.shape != ct.hu.shape:
        raise ValueError("Structure mask shape does not match CT.")
    if not np.any(selected):
        raise ValueError("Structure mask contains no voxels.")

    orientation = np.asarray(
        ct.geometry.image_orientation_patient,
        dtype=np.float64,
    )
    column_direction = orientation[:3]
    row_direction = orientation[3:]
    row_spacing, column_spacing = ct.geometry.pixel_spacing_mm

    chunks: list[np.ndarray] = []
    for frame in range(selected.shape[0]):
        rows, columns = np.nonzero(selected[frame])
        if rows.size == 0:
            continue

        origin = ct.image_positions_patient_mm[frame]
        points = (
            origin[None, :]
            + columns[:, None] * column_spacing * column_direction[None, :]
            + rows[:, None] * row_spacing * row_direction[None, :]
        )
        chunks.append(points)

    return np.concatenate(chunks, axis=0)


def sample_structure_dose_gy(
    ct: CtSeries,
    structure_mask: np.ndarray,
    sampler: PatientDoseSamplerGy,
) -> StructureDoseSamples:
    """Sample dose at CT voxel centres inside a rasterized structure."""

    if sampler.frame_of_reference_uid != ct.geometry.frame_of_reference_uid:
        raise ValueError(
            "Dose sampler and CT use different FrameOfReferenceUID values."
        )

    points = _ct_mask_patient_points(ct, structure_mask)
    dose = np.asarray(
        sampler.sample_patient_points_gy(points),
        dtype=np.float64,
    )
    if dose.shape != (points.shape[0],):
        raise ValueError("Dose sampler returned unexpected structure-dose shape.")
    if not np.all(np.isfinite(dose)):
        raise ValueError("Structure dose contains non-finite values.")

    return StructureDoseSamples(
        dose_gy=dose,
        voxel_volume_cm3=ct.geometry.voxel_volume_mm3 / 1000.0,
    )


def dvh_metrics(samples: StructureDoseSamples) -> DvhMetrics:
    dose = np.asarray(samples.dose_gy, dtype=np.float64)
    if dose.ndim != 1 or dose.size == 0:
        raise ValueError("DVH dose samples must be a non-empty 1-D array.")
    if not np.all(np.isfinite(dose)):
        raise ValueError("DVH dose samples must be finite.")

    return DvhMetrics(
        volume_cm3=samples.volume_cm3,
        d98_gy=float(np.percentile(dose, 2)),
        d95_gy=float(np.percentile(dose, 5)),
        d50_gy=float(np.percentile(dose, 50)),
        d2_gy=float(np.percentile(dose, 98)),
        mean_gy=float(np.mean(dose)),
        min_gy=float(np.min(dose)),
        max_gy=float(np.max(dose)),
    )
