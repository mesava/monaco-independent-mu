from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from indep_mu.dicom.ct import CtSeries

from .density import red_to_mass_density
from .hu_red import CtCalibration, DRT120KV
from .materials import (
    MaterialBlendArrays,
    PATIENT_LUT,
    material_blend_arrays_for_density,
)


@dataclass(frozen=True)
class CropBounds:
    z_start: int
    z_stop: int
    y_start: int
    y_stop: int
    x_start: int
    x_stop: int

    @property
    def slices(self) -> tuple[slice, slice, slice]:
        return (
            slice(self.z_start, self.z_stop),
            slice(self.y_start, self.y_stop),
            slice(self.x_start, self.x_stop),
        )


@dataclass(frozen=True)
class PatientModel:
    """Voxel-wise Monaco-matched patient model on a cropped CT grid."""

    hu: np.ndarray
    red: np.ndarray
    mass_density_g_cm3: np.ndarray
    patient_mask: np.ndarray
    out_of_range_code: np.ndarray
    material_blends: MaterialBlendArrays
    crop_bounds: CropBounds
    voxel_volume_mm3: float

    @property
    def shape(self) -> tuple[int, int, int]:
        return tuple(int(x) for x in self.hu.shape)

    @property
    def patient_voxel_count(self) -> int:
        return int(np.count_nonzero(self.patient_mask))

    @property
    def patient_volume_cm3(self) -> float:
        return self.patient_voxel_count * self.voxel_volume_mm3 / 1000.0

    @property
    def low_hu_count(self) -> int:
        return int(np.count_nonzero(self.out_of_range_code == -1))

    @property
    def high_hu_count(self) -> int:
        return int(np.count_nonzero(self.out_of_range_code == 1))


def crop_bounds_from_mask(
    mask: np.ndarray,
    *,
    margin_zyx: tuple[int, int, int] = (0, 0, 0),
) -> CropBounds:
    data = np.asarray(mask, dtype=bool)
    if data.ndim != 3:
        raise ValueError("Patient mask must be a 3-D array.")
    if not np.any(data):
        raise ValueError("Patient mask contains no voxels.")
    if any(margin < 0 for margin in margin_zyx):
        raise ValueError("Crop margins must be non-negative.")

    coordinates = np.argwhere(data)
    minimum = coordinates.min(axis=0)
    maximum = coordinates.max(axis=0) + 1

    starts = np.maximum(minimum - np.asarray(margin_zyx), 0)
    stops = np.minimum(maximum + np.asarray(margin_zyx), data.shape)

    return CropBounds(
        z_start=int(starts[0]),
        z_stop=int(stops[0]),
        y_start=int(starts[1]),
        y_stop=int(stops[1]),
        x_start=int(starts[2]),
        x_stop=int(stops[2]),
    )


def build_patient_model(
    ct: CtSeries,
    patient_mask: np.ndarray,
    *,
    calibration: CtCalibration = DRT120KV,
    out_of_range: str = "clip",
    margin_zyx: tuple[int, int, int] = (0, 0, 0),
) -> PatientModel:
    """Build HU -> RED -> density -> material/material-mix arrays.

    Out-of-range HU is evaluated only inside the Patient ROI. Outside the
    Patient ROI the voxel is explicitly treated as dry air using the minimum
    RED from the scanner calibration.

    out_of_range:
      - "raise": fail if Patient voxels leave the scanner calibration range.
      - "clip": reproduce Monaco endpoint RED behaviour while recording every
        clipped Patient voxel in out_of_range_code.
    """

    mask = np.asarray(patient_mask, dtype=bool)
    if mask.shape != ct.hu.shape:
        raise ValueError(
            f"Patient mask shape {mask.shape} does not match CT shape {ct.hu.shape}."
        )
    if out_of_range not in {"raise", "clip"}:
        raise ValueError("out_of_range must be 'raise' or 'clip'.")

    bounds = crop_bounds_from_mask(mask, margin_zyx=margin_zyx)
    slc = bounds.slices

    hu = np.asarray(ct.hu[slc], dtype=np.float32).copy()
    body = np.asarray(mask[slc], dtype=bool).copy()

    below = body & (hu < calibration.hu_min)
    above = body & (hu > calibration.hu_max)

    if out_of_range == "raise" and (np.any(below) or np.any(above)):
        raise ValueError(
            "Patient ROI contains HU values outside CT calibration range: "
            f"below={np.count_nonzero(below)}, above={np.count_nonzero(above)}."
        )

    # Monaco-matched conversion uses endpoint clipping. We always do that
    # internally after the safety check so voxels outside the external Patient
    # contour cannot fail the patient calculation because of scanner padding.
    red = calibration.to_red(hu, out_of_range="clip").astype(np.float32)

    # Voxels outside Patient are not inferred from arbitrary CT background;
    # they are explicitly dry air. PPS/couch is a separate geometry layer.
    red[~body] = np.float32(calibration.red[0])

    density = red_to_mass_density(red).astype(np.float32)

    blend = material_blend_arrays_for_density(density, PATIENT_LUT)

    out_code = np.zeros(hu.shape, dtype=np.int8)
    out_code[below] = -1
    out_code[above] = 1

    return PatientModel(
        hu=hu,
        red=red,
        mass_density_g_cm3=density,
        patient_mask=body,
        out_of_range_code=out_code,
        material_blends=blend,
        crop_bounds=bounds,
        voxel_volume_mm3=ct.geometry.voxel_volume_mm3,
    )
