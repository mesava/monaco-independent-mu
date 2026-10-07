from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol, runtime_checkable

import numpy as np

from .dosegrid import RectilinearDoseGrid, trilinear_sample


@dataclass(frozen=True)
class SphericalRoi:
    name: str
    center_cm: tuple[float, float, float]
    radius_cm: float = 0.25

    def __post_init__(self) -> None:
        centre = np.asarray(self.center_cm, dtype=np.float64)
        if centre.shape != (3,) or not np.all(np.isfinite(centre)):
            raise ValueError("ROI center must contain three finite coordinates.")
        if self.radius_cm <= 0 or not np.isfinite(self.radius_cm):
            raise ValueError("ROI radius must be finite and positive.")


@dataclass(frozen=True)
class RoiDoseStatistic:
    roi: SphericalRoi
    mean_dose: float
    sample_count: int


def _sphere_sample_points(
    roi: SphericalRoi,
    *,
    samples_per_axis: int,
) -> np.ndarray:
    if samples_per_axis < 3 or samples_per_axis % 2 == 0:
        raise ValueError("samples_per_axis must be an odd integer >= 3.")

    offsets = np.linspace(
        -roi.radius_cm,
        roi.radius_cm,
        samples_per_axis,
        dtype=np.float64,
    )
    xx, yy, zz = np.meshgrid(offsets, offsets, offsets, indexing="xy")
    points = np.column_stack([xx.ravel(), yy.ravel(), zz.ravel()])
    inside = np.sum(points**2, axis=1) <= roi.radius_cm**2 + 1e-15
    return points[inside] + np.asarray(roi.center_cm, dtype=np.float64)


def spherical_roi_mean(
    grid: RectilinearDoseGrid,
    roi: SphericalRoi,
    *,
    samples_per_axis: int = 11,
) -> RoiDoseStatistic:
    """Approximate the volume mean with symmetric supersampling.

    Sampling is performed in continuous coordinates followed by trilinear
    interpolation. This avoids defining a 2.5-mm-radius ROI as a fragile set of
    whichever raw dose voxels happen to contain their centres.
    """

    points = _sphere_sample_points(roi, samples_per_axis=samples_per_axis)
    values = trilinear_sample(grid, points)
    return RoiDoseStatistic(
        roi=roi,
        mean_dose=float(np.mean(values)),
        sample_count=int(values.size),
    )


@runtime_checkable
class PatientDoseSamplerGy(Protocol):
    """Dose source that samples Gy at DICOM patient coordinates in mm."""

    def sample_patient_points_gy(self, points_patient_mm: np.ndarray) -> np.ndarray:
        ...


@dataclass(frozen=True)
class PatientSphericalRoi:
    """Spherical ROI whose centre is expressed in DICOM patient coordinates."""

    name: str
    center_patient_mm: tuple[float, float, float]
    radius_mm: float = 2.5

    def __post_init__(self) -> None:
        centre = np.asarray(self.center_patient_mm, dtype=np.float64)
        if centre.shape != (3,) or not np.all(np.isfinite(centre)):
            raise ValueError(
                "Patient ROI center must contain three finite coordinates."
            )
        if self.radius_mm <= 0 or not np.isfinite(self.radius_mm):
            raise ValueError("Patient ROI radius must be finite and positive.")


@dataclass(frozen=True)
class PatientRoiDoseStatistic:
    roi: PatientSphericalRoi
    mean_dose_gy: float
    sample_count: int


def _patient_sphere_sample_points(
    roi: PatientSphericalRoi,
    *,
    samples_per_axis: int,
) -> np.ndarray:
    if samples_per_axis < 3 or samples_per_axis % 2 == 0:
        raise ValueError("samples_per_axis must be an odd integer >= 3.")

    offsets = np.linspace(
        -roi.radius_mm,
        roi.radius_mm,
        samples_per_axis,
        dtype=np.float64,
    )
    xx, yy, zz = np.meshgrid(offsets, offsets, offsets, indexing="xy")
    points = np.column_stack([xx.ravel(), yy.ravel(), zz.ravel()])
    inside = np.sum(points**2, axis=1) <= roi.radius_mm**2 + 1e-12
    return points[inside] + np.asarray(roi.center_patient_mm, dtype=np.float64)


def spherical_patient_roi_mean_gy(
    sampler: PatientDoseSamplerGy,
    roi: PatientSphericalRoi,
    *,
    samples_per_axis: int = 11,
) -> PatientRoiDoseStatistic:
    """Mean dose in a patient-coordinate sphere using continuous sampling.

    RTDOSE and MC adapters can implement the same sampler protocol. This keeps
    ROI coordinates identical for TPS and independent dose distributions and
    prevents accidental comparison of different coordinate systems.
    """

    points = _patient_sphere_sample_points(
        roi,
        samples_per_axis=samples_per_axis,
    )
    values = np.asarray(
        sampler.sample_patient_points_gy(points),
        dtype=np.float64,
    )
    if values.shape != (points.shape[0],):
        raise ValueError(
            "Patient dose sampler returned an unexpected result shape."
        )
    if not np.all(np.isfinite(values)):
        raise ValueError("Patient dose sampler returned non-finite values.")

    return PatientRoiDoseStatistic(
        roi=roi,
        mean_dose_gy=float(np.mean(values)),
        sample_count=int(values.size),
    )
