from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from indep_mu.dicom.ct import CtSeries


@dataclass(frozen=True)
class PatientHuSummary:
    voxel_count: int
    volume_cm3: float
    hu_min: float
    hu_max: float
    hu_mean: float
    below_calibration_count: int
    above_calibration_count: int

    @property
    def below_calibration_fraction(self) -> float:
        return self.below_calibration_count / self.voxel_count

    @property
    def above_calibration_fraction(self) -> float:
        return self.above_calibration_count / self.voxel_count


def summarize_patient_hu(
    ct: CtSeries,
    patient_mask: np.ndarray,
    *,
    calibration_hu_min: float,
    calibration_hu_max: float,
) -> PatientHuSummary:
    mask = np.asarray(patient_mask, dtype=bool)
    if mask.shape != ct.hu.shape:
        raise ValueError(
            f"Patient mask shape {mask.shape} does not match CT shape {ct.hu.shape}."
        )

    values = ct.hu[mask]
    if values.size == 0:
        raise ValueError("Patient mask contains no voxels.")

    voxel_volume_cm3 = ct.geometry.voxel_volume_mm3 / 1000.0

    return PatientHuSummary(
        voxel_count=int(values.size),
        volume_cm3=float(values.size * voxel_volume_cm3),
        hu_min=float(np.min(values)),
        hu_max=float(np.max(values)),
        hu_mean=float(np.mean(values)),
        below_calibration_count=int(np.count_nonzero(values < calibration_hu_min)),
        above_calibration_count=int(np.count_nonzero(values > calibration_hu_max)),
    )
