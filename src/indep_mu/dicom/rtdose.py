from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pydicom


@dataclass(frozen=True)
class DoseGeometry:
    rows: int
    columns: int
    frames: int
    pixel_spacing_mm: tuple[float, float]
    image_orientation_patient: tuple[float, float, float, float, float, float]
    image_position_patient_mm: tuple[float, float, float]
    frame_offsets_mm: np.ndarray
    frame_positions_patient_mm: np.ndarray
    frame_offset_mode: str

    @property
    def column_direction(self) -> np.ndarray:
        return np.asarray(self.image_orientation_patient[:3], dtype=np.float64)

    @property
    def row_direction(self) -> np.ndarray:
        return np.asarray(self.image_orientation_patient[3:], dtype=np.float64)

    @property
    def normal_direction(self) -> np.ndarray:
        normal = np.cross(self.column_direction, self.row_direction)
        return normal / np.linalg.norm(normal)

    def voxel_center_patient_mm(
        self,
        frame: int,
        row: int,
        column: int,
    ) -> np.ndarray:
        if not 0 <= frame < self.frames:
            raise IndexError("frame out of range")
        if not 0 <= row < self.rows:
            raise IndexError("row out of range")
        if not 0 <= column < self.columns:
            raise IndexError("column out of range")

        row_spacing, column_spacing = self.pixel_spacing_mm
        return (
            self.frame_positions_patient_mm[frame]
            + column * column_spacing * self.column_direction
            + row * row_spacing * self.row_direction
        )


@dataclass(frozen=True)
class RtDose:
    dose: np.ndarray
    geometry: DoseGeometry
    dose_units: str
    dose_type: str
    dose_summation_type: str
    frame_of_reference_uid: str
    sop_instance_uid: str
    referenced_rtplan_uids: tuple[str, ...]
    referenced_beam_numbers: tuple[int, ...]

    @property
    def shape(self) -> tuple[int, int, int]:
        return tuple(int(value) for value in self.dose.shape)


def _float_tuple(value, length: int, name: str) -> tuple[float, ...]:
    result = tuple(float(item) for item in value)
    if len(result) != length:
        raise ValueError(f"{name} must contain {length} values.")
    return result


def _frame_positions(
    image_position: np.ndarray,
    orientation: tuple[float, float, float, float, float, float],
    offsets: np.ndarray,
    *,
    tolerance_mm: float,
) -> tuple[np.ndarray, str]:
    column_direction = np.asarray(orientation[:3], dtype=np.float64)
    row_direction = np.asarray(orientation[3:], dtype=np.float64)
    normal = np.cross(column_direction, row_direction)
    norm = float(np.linalg.norm(normal))
    if norm == 0:
        raise ValueError("Invalid ImageOrientationPatient: zero dose-grid normal.")
    normal /= norm

    if np.isclose(offsets[0], 0.0, rtol=0.0, atol=tolerance_mm):
        return (
            image_position[None, :] + offsets[:, None] * normal[None, :],
            "RELATIVE",
        )

    axial = np.allclose(
        np.asarray(orientation, dtype=np.float64),
        np.asarray([1.0, 0.0, 0.0, 0.0, 1.0, 0.0]),
        rtol=0.0,
        atol=1e-7,
    )
    if axial and np.isclose(
        offsets[0],
        image_position[2],
        rtol=0.0,
        atol=tolerance_mm,
    ):
        positions = np.repeat(image_position[None, :], offsets.size, axis=0)
        positions[:, 2] = offsets
        return positions, "ABSOLUTE_PATIENT_Z_LEGACY"

    raise ValueError(
        "GridFrameOffsetVector is neither DICOM relative option (first value 0) "
        "nor legacy absolute patient-z option for axial orientation."
    )


def _referenced_plan_and_beams(
    dataset: pydicom.dataset.Dataset,
) -> tuple[tuple[str, ...], tuple[int, ...]]:
    plan_uids: list[str] = []
    beam_numbers: set[int] = set()

    for plan_ref in getattr(dataset, "ReferencedRTPlanSequence", []) or []:
        uid = getattr(plan_ref, "ReferencedSOPInstanceUID", None)
        if uid is not None:
            plan_uids.append(str(uid))

        for beam_ref in getattr(plan_ref, "ReferencedBeamSequence", []) or []:
            if hasattr(beam_ref, "ReferencedBeamNumber"):
                beam_numbers.add(int(beam_ref.ReferencedBeamNumber))

        for fraction_ref in (
            getattr(plan_ref, "ReferencedFractionGroupSequence", []) or []
        ):
            for beam_ref in (
                getattr(fraction_ref, "ReferencedBeamSequence", []) or []
            ):
                if hasattr(beam_ref, "ReferencedBeamNumber"):
                    beam_numbers.add(int(beam_ref.ReferencedBeamNumber))

    for beam_ref in getattr(dataset, "ReferencedBeamSequence", []) or []:
        if hasattr(beam_ref, "ReferencedBeamNumber"):
            beam_numbers.add(int(beam_ref.ReferencedBeamNumber))

    return tuple(plan_uids), tuple(sorted(beam_numbers))


def load_rtdose(
    path: str | Path,
    *,
    geometry_tolerance_mm: float = 1e-4,
) -> RtDose:
    """Load RT Dose and construct dose-plane coordinates in patient space.

    GridFrameOffsetVector follows DICOM PS3.3 C.8.8.3.2:
    - first element 0 -> offsets relative to ImagePositionPatient along the
      cross-product of image row/column directions;
    - legacy absolute patient-z coordinates are accepted only for unrotated
      transverse orientation.
    """

    dataset = pydicom.dcmread(path, force=False)
    if getattr(dataset, "Modality", None) != "RTDOSE":
        raise ValueError("Input file is not RT Dose.")

    required = (
        "Rows",
        "Columns",
        "PixelSpacing",
        "ImageOrientationPatient",
        "ImagePositionPatient",
        "DoseGridScaling",
        "DoseUnits",
        "DoseType",
        "DoseSummationType",
        "FrameOfReferenceUID",
        "SOPInstanceUID",
    )
    missing = [name for name in required if not hasattr(dataset, name)]
    if missing:
        raise ValueError(
            "RTDOSE is missing required attributes: " + ", ".join(missing)
        )

    rows = int(dataset.Rows)
    columns = int(dataset.Columns)
    frames = int(getattr(dataset, "NumberOfFrames", 1))
    pixel_spacing = _float_tuple(dataset.PixelSpacing, 2, "PixelSpacing")
    orientation = _float_tuple(
        dataset.ImageOrientationPatient,
        6,
        "ImageOrientationPatient",
    )
    image_position = np.asarray(
        _float_tuple(dataset.ImagePositionPatient, 3, "ImagePositionPatient"),
        dtype=np.float64,
    )

    if frames > 1:
        if not hasattr(dataset, "GridFrameOffsetVector"):
            raise ValueError(
                "Multi-frame RTDOSE has no GridFrameOffsetVector."
            )
        offsets = np.asarray(
            [float(value) for value in dataset.GridFrameOffsetVector],
            dtype=np.float64,
        )
        if offsets.size != frames:
            raise ValueError(
                f"GridFrameOffsetVector has {offsets.size} values, "
                f"NumberOfFrames={frames}."
            )
    else:
        offsets = np.asarray([0.0], dtype=np.float64)

    if offsets.size > 1:
        differences = np.diff(offsets)
        if not (np.all(differences > 0) or np.all(differences < 0)):
            raise ValueError("GridFrameOffsetVector must be strictly monotonic.")

    frame_positions, offset_mode = _frame_positions(
        image_position,
        orientation,
        offsets,
        tolerance_mm=geometry_tolerance_mm,
    )

    pixels = np.asarray(dataset.pixel_array)
    if frames == 1:
        if pixels.shape != (rows, columns):
            raise ValueError(
                f"RTDOSE pixel shape {pixels.shape} != {(rows, columns)}."
            )
        pixels = pixels[np.newaxis, :, :]
    elif pixels.shape != (frames, rows, columns):
        raise ValueError(
            f"RTDOSE pixel shape {pixels.shape} != "
            f"{(frames, rows, columns)}."
        )

    dose = (
        pixels.astype(np.float64)
        * float(dataset.DoseGridScaling)
    ).astype(np.float32)

    plan_uids, beam_numbers = _referenced_plan_and_beams(dataset)

    geometry = DoseGeometry(
        rows=rows,
        columns=columns,
        frames=frames,
        pixel_spacing_mm=(float(pixel_spacing[0]), float(pixel_spacing[1])),
        image_orientation_patient=tuple(float(v) for v in orientation),
        image_position_patient_mm=tuple(float(v) for v in image_position),
        frame_offsets_mm=offsets,
        frame_positions_patient_mm=frame_positions,
        frame_offset_mode=offset_mode,
    )

    return RtDose(
        dose=dose,
        geometry=geometry,
        dose_units=str(dataset.DoseUnits),
        dose_type=str(dataset.DoseType),
        dose_summation_type=str(dataset.DoseSummationType),
        frame_of_reference_uid=str(dataset.FrameOfReferenceUID),
        sop_instance_uid=str(dataset.SOPInstanceUID),
        referenced_rtplan_uids=plan_uids,
        referenced_beam_numbers=beam_numbers,
    )
