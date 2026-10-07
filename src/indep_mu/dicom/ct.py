from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pydicom


@dataclass(frozen=True)
class CtGeometry:
    """Geometry of a validated CT image stack in the DICOM patient frame."""

    rows: int
    columns: int
    pixel_spacing_mm: tuple[float, float]
    slice_spacing_mm: float
    image_orientation_patient: tuple[float, float, float, float, float, float]
    first_image_position_patient_mm: tuple[float, float, float]
    last_image_position_patient_mm: tuple[float, float, float]
    patient_position: str | None
    frame_of_reference_uid: str
    series_instance_uid: str
    kvp: float | None

    @property
    def voxel_volume_mm3(self) -> float:
        return (
            self.pixel_spacing_mm[0]
            * self.pixel_spacing_mm[1]
            * self.slice_spacing_mm
        )


@dataclass(frozen=True)
class CtSeries:
    """Validated CT volume and geometry.

    The HU array is ordered along the slice-normal direction and has shape
    (slices, rows, columns).
    """

    hu: np.ndarray
    geometry: CtGeometry
    slice_positions_mm: np.ndarray
    image_positions_patient_mm: np.ndarray
    sop_instance_uids: tuple[str, ...]

    def calibration_range_summary(
        self,
        hu_min: float,
        hu_max: float,
        *,
        mask: np.ndarray | None = None,
    ) -> dict[str, int | float]:
        values = self.hu if mask is None else self.hu[np.asarray(mask, dtype=bool)]
        if values.size == 0:
            raise ValueError("Calibration summary mask contains no voxels.")

        below = values < hu_min
        above = values > hu_max
        return {
            "voxel_count": int(values.size),
            "hu_observed_min": float(np.min(values)),
            "hu_observed_max": float(np.max(values)),
            "hu_observed_mean": float(np.mean(values)),
            "below_calibration_count": int(np.count_nonzero(below)),
            "above_calibration_count": int(np.count_nonzero(above)),
        }


def _as_float_tuple(value: Iterable[float], expected_length: int) -> tuple[float, ...]:
    result = tuple(float(x) for x in value)
    if len(result) != expected_length:
        raise ValueError(
            f"Expected {expected_length} values, received {len(result)}."
        )
    return result


def _slice_normal(
    orientation: tuple[float, float, float, float, float, float],
) -> np.ndarray:
    column_index_direction = np.asarray(orientation[:3], dtype=np.float64)
    row_index_direction = np.asarray(orientation[3:], dtype=np.float64)
    normal = np.cross(column_index_direction, row_index_direction)
    norm = float(np.linalg.norm(normal))
    if norm == 0:
        raise ValueError("Invalid ImageOrientationPatient: zero slice normal.")
    return normal / norm


def _slice_coordinate(ds: pydicom.dataset.Dataset, normal: np.ndarray) -> float:
    position = np.asarray(
        _as_float_tuple(ds.ImagePositionPatient, 3), dtype=np.float64
    )
    return float(np.dot(position, normal))


def _close_tuple(
    a: tuple[float, ...],
    b: tuple[float, ...],
    *,
    atol: float,
) -> bool:
    return bool(np.allclose(a, b, rtol=0.0, atol=atol))


def load_ct_series(
    directory: str | Path,
    *,
    series_instance_uid: str | None = None,
    geometry_tolerance_mm: float = 1e-4,
    spacing_tolerance_mm: float = 1e-3,
) -> CtSeries:
    """Load and validate one CT DICOM series from a directory.

    When series_instance_uid is supplied, unrelated CT series under the same
    export directory are ignored. Without it, the loader deliberately fails on
    mixed SeriesInstanceUID values rather than guessing which CT is clinically
    referenced.

    The selected series is also checked for mixed FrameOfReferenceUID,
    inconsistent matrix/pixel spacing/orientation, duplicate slice locations,
    and non-uniform slice spacing.

    Patient names/IDs are intentionally not retained by this data model.
    """

    root = Path(directory)
    if not root.exists():
        raise FileNotFoundError(root)

    datasets: list[pydicom.dataset.Dataset] = []
    for path in sorted(p for p in root.rglob("*") if p.is_file()):
        try:
            ds = pydicom.dcmread(path, stop_before_pixels=True, force=False)
        except Exception:
            continue
        if getattr(ds, "Modality", None) == "CT":
            if (
                series_instance_uid is not None
                and str(getattr(ds, "SeriesInstanceUID", ""))
                != str(series_instance_uid)
            ):
                continue
            ds.filename = str(path)
            datasets.append(ds)

    if not datasets:
        if series_instance_uid is None:
            raise ValueError(f"No CT DICOM instances found under {root}.")
        raise ValueError(
            "No CT DICOM instances found for requested SeriesInstanceUID "
            f"{series_instance_uid!r} under {root}."
        )

    reference = datasets[0]

    required = (
        "Rows",
        "Columns",
        "PixelSpacing",
        "ImagePositionPatient",
        "ImageOrientationPatient",
        "FrameOfReferenceUID",
        "SeriesInstanceUID",
        "SOPInstanceUID",
    )
    for ds in datasets:
        missing = [name for name in required if not hasattr(ds, name)]
        if missing:
            raise ValueError(
                f"CT instance {getattr(ds, 'filename', '<unknown>')} "
                f"is missing required attributes: {', '.join(missing)}"
            )

    rows = int(reference.Rows)
    columns = int(reference.Columns)
    pixel_spacing = _as_float_tuple(reference.PixelSpacing, 2)
    orientation = _as_float_tuple(reference.ImageOrientationPatient, 6)
    frame_uid = str(reference.FrameOfReferenceUID)
    series_uid = str(reference.SeriesInstanceUID)
    normal = _slice_normal(orientation)

    for ds in datasets:
        if int(ds.Rows) != rows or int(ds.Columns) != columns:
            raise ValueError("CT matrix size is inconsistent across slices.")
        if str(ds.FrameOfReferenceUID) != frame_uid:
            raise ValueError("Multiple FrameOfReferenceUID values in CT input.")
        if str(ds.SeriesInstanceUID) != series_uid:
            raise ValueError("Multiple CT SeriesInstanceUID values in input.")
        if not _close_tuple(
            _as_float_tuple(ds.PixelSpacing, 2),
            pixel_spacing,
            atol=geometry_tolerance_mm,
        ):
            raise ValueError("PixelSpacing is inconsistent across CT slices.")
        if not _close_tuple(
            _as_float_tuple(ds.ImageOrientationPatient, 6),
            orientation,
            atol=1e-6,
        ):
            raise ValueError(
                "ImageOrientationPatient is inconsistent across CT slices."
            )

    ordered = sorted(datasets, key=lambda ds: _slice_coordinate(ds, normal))
    positions = np.asarray(
        [_slice_coordinate(ds, normal) for ds in ordered], dtype=np.float64
    )
    image_positions = np.asarray(
        [
            _as_float_tuple(ds.ImagePositionPatient, 3)
            for ds in ordered
        ],
        dtype=np.float64,
    )

    if len(positions) < 2:
        raise ValueError("At least two CT slices are required for a 3D volume.")

    differences = np.diff(positions)
    if np.any(differences <= 0):
        raise ValueError("Duplicate or non-monotonic CT slice positions.")

    slice_spacing = float(np.median(differences))
    if not np.allclose(
        differences,
        slice_spacing,
        rtol=0.0,
        atol=spacing_tolerance_mm,
    ):
        observed = (float(np.min(differences)), float(np.max(differences)))
        raise ValueError(
            "Non-uniform CT slice spacing: "
            f"median={slice_spacing:.6f} mm, observed={observed}."
        )

    hu_slices: list[np.ndarray] = []
    sop_uids: list[str] = []

    for header in ordered:
        ds = pydicom.dcmread(header.filename, force=False)
        pixels = np.asarray(ds.pixel_array, dtype=np.float32)
        if pixels.shape != (rows, columns):
            raise ValueError(
                f"Unexpected decoded pixel matrix {pixels.shape}; "
                f"expected {(rows, columns)}."
            )

        slope = float(getattr(ds, "RescaleSlope", 1.0))
        intercept = float(getattr(ds, "RescaleIntercept", 0.0))
        hu_slices.append(pixels * slope + intercept)
        sop_uids.append(str(ds.SOPInstanceUID))

    hu = np.stack(hu_slices, axis=0)

    first_position = tuple(float(x) for x in image_positions[0])
    last_position = tuple(float(x) for x in image_positions[-1])

    patient_position = getattr(reference, "PatientPosition", None)
    kvp_value = getattr(reference, "KVP", None)

    geometry = CtGeometry(
        rows=rows,
        columns=columns,
        pixel_spacing_mm=(float(pixel_spacing[0]), float(pixel_spacing[1])),
        slice_spacing_mm=slice_spacing,
        image_orientation_patient=tuple(float(x) for x in orientation),
        first_image_position_patient_mm=first_position,
        last_image_position_patient_mm=last_position,
        patient_position=str(patient_position) if patient_position else None,
        frame_of_reference_uid=frame_uid,
        series_instance_uid=series_uid,
        kvp=float(kvp_value) if kvp_value is not None else None,
    )

    return CtSeries(
        hu=hu,
        geometry=geometry,
        slice_positions_mm=positions,
        image_positions_patient_mm=image_positions,
        sop_instance_uids=tuple(sop_uids),
    )
