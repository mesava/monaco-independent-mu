from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pydicom

from .ct import CtSeries


@dataclass(frozen=True)
class RoiDefinition:
    roi_number: int
    roi_name: str
    frame_of_reference_uid: str | None


@dataclass(frozen=True)
class RoiMask:
    """Rasterized RTSTRUCT ROI on the native CT voxel grid."""

    roi_number: int
    roi_name: str
    mask: np.ndarray
    contour_count: int
    referenced_slice_count: int

    @property
    def voxel_count(self) -> int:
        return int(np.count_nonzero(self.mask))


def list_rtstruct_rois(rtstruct_path: str | Path) -> tuple[RoiDefinition, ...]:
    """List RTSTRUCT ROI definitions without rasterizing contours."""

    rtstruct = pydicom.dcmread(rtstruct_path, stop_before_pixels=True, force=False)
    if getattr(rtstruct, "Modality", None) != "RTSTRUCT":
        raise ValueError("Input file is not an RT Structure Set.")

    sequence = getattr(rtstruct, "StructureSetROISequence", None)
    if not sequence:
        return ()

    result: list[RoiDefinition] = []
    seen_numbers: set[int] = set()
    for item in sequence:
        number = int(item.ROINumber)
        if number in seen_numbers:
            raise ValueError(f"Duplicate ROINumber in RTSTRUCT: {number}.")
        seen_numbers.add(number)
        frame = getattr(item, "ReferencedFrameOfReferenceUID", None)
        result.append(
            RoiDefinition(
                roi_number=number,
                roi_name=str(item.ROIName),
                frame_of_reference_uid=str(frame) if frame is not None else None,
            )
        )

    return tuple(result)


def _roi_number_by_name(rtstruct: pydicom.dataset.Dataset, roi_name: str) -> int:
    matches = [
        int(item.ROINumber)
        for item in rtstruct.StructureSetROISequence
        if str(item.ROIName).casefold() == roi_name.casefold()
    ]
    if not matches:
        raise ValueError(f"ROI not found in RTSTRUCT: {roi_name!r}")
    if len(matches) > 1:
        raise ValueError(f"ROI name is not unique in RTSTRUCT: {roi_name!r}")
    return matches[0]


def _polygon_mask(
    rows: int,
    columns: int,
    row_coordinates: np.ndarray,
    column_coordinates: np.ndarray,
) -> np.ndarray:
    """Rasterize a polygon using an even-odd scanline rule at voxel centres."""

    result = np.zeros((rows, columns), dtype=bool)
    if row_coordinates.size < 3:
        return result

    r = np.asarray(row_coordinates, dtype=np.float64)
    c = np.asarray(column_coordinates, dtype=np.float64)

    r_next = np.roll(r, -1)
    c_next = np.roll(c, -1)

    row_min = max(0, int(np.ceil(np.min(r))))
    row_max = min(rows - 1, int(np.floor(np.max(r))))

    for row_index in range(row_min, row_max + 1):
        y = float(row_index)
        crosses = ((r <= y) & (r_next > y)) | ((r_next <= y) & (r > y))
        if not np.any(crosses):
            continue

        r0 = r[crosses]
        r1 = r_next[crosses]
        c0 = c[crosses]
        c1 = c_next[crosses]
        intersections = c0 + (y - r0) * (c1 - c0) / (r1 - r0)
        intersections.sort()

        if intersections.size % 2 != 0:
            raise ValueError(
                "Polygon scanline produced an odd number of intersections; "
                "contour may be self-intersecting."
            )

        for left, right in intersections.reshape(-1, 2):
            first = max(0, int(np.ceil(min(left, right))))
            last = min(columns - 1, int(np.floor(max(left, right))))
            if first <= last:
                result[row_index, first : last + 1] = True

    return result


def build_roi_mask(
    ct: CtSeries,
    rtstruct_path: str | Path,
    roi_name: str,
    *,
    contour_plane_tolerance_mm: float = 0.25,
) -> RoiMask:
    """Rasterize one RTSTRUCT ROI onto the native CT grid.

    Contours are attached to CT slices through ContourImageSequence /
    ReferencedSOPInstanceUID. Silent nearest-slice matching is deliberately not
    used: a missing referenced CT instance is a hard error.

    CLOSED_PLANAR contours are unioned. CLOSEDPLANAR_XOR contours are combined
    with XOR, as required by their geometric type.
    """

    rtstruct = pydicom.dcmread(rtstruct_path, force=False)
    if getattr(rtstruct, "Modality", None) != "RTSTRUCT":
        raise ValueError("Input file is not an RT Structure Set.")

    roi_number = _roi_number_by_name(rtstruct, roi_name)

    roi_definition = next(
        item
        for item in rtstruct.StructureSetROISequence
        if int(item.ROINumber) == roi_number
    )
    referenced_frame = getattr(roi_definition, "ReferencedFrameOfReferenceUID", None)
    if (
        referenced_frame is not None
        and str(referenced_frame) != ct.geometry.frame_of_reference_uid
    ):
        raise ValueError(
            "RTSTRUCT ROI and CT have different FrameOfReferenceUID values."
        )

    roi_contour = next(
        (
            item
            for item in rtstruct.ROIContourSequence
            if int(item.ReferencedROINumber) == roi_number
        ),
        None,
    )
    if roi_contour is None:
        raise ValueError(f"ROI {roi_name!r} has no ROIContourSequence entry.")

    uid_to_slice = {
        uid: index for index, uid in enumerate(ct.sop_instance_uids)
    }

    orientation = np.asarray(
        ct.geometry.image_orientation_patient, dtype=np.float64
    )
    column_index_direction = orientation[:3]
    row_index_direction = orientation[3:]
    slice_normal = np.cross(column_index_direction, row_index_direction)
    slice_normal /= np.linalg.norm(slice_normal)

    row_spacing, column_spacing = ct.geometry.pixel_spacing_mm
    mask = np.zeros_like(ct.hu, dtype=bool)

    contour_sequence = getattr(roi_contour, "ContourSequence", [])
    contour_count = 0
    referenced_slices: set[int] = set()
    geometric_types: set[str] = set()

    for contour in contour_sequence:
        geometric_type = str(getattr(contour, "ContourGeometricType", ""))
        if geometric_type not in {"CLOSED_PLANAR", "CLOSEDPLANAR_XOR"}:
            raise ValueError(
                f"Unsupported contour geometric type for {roi_name!r}: "
                f"{geometric_type!r}"
            )
        geometric_types.add(geometric_type)

        image_sequence = getattr(contour, "ContourImageSequence", None)
        if not image_sequence or len(image_sequence) != 1:
            raise ValueError(
                f"Contour in ROI {roi_name!r} does not reference exactly one CT slice."
            )

        referenced_uid = str(image_sequence[0].ReferencedSOPInstanceUID)
        if referenced_uid not in uid_to_slice:
            raise ValueError(
                f"RTSTRUCT contour references CT SOPInstanceUID not present "
                f"in the loaded CT series: {referenced_uid}"
            )

        slice_index = uid_to_slice[referenced_uid]
        image_position = ct.image_positions_patient_mm[slice_index]

        data = np.asarray(contour.ContourData, dtype=np.float64)
        if data.size % 3 != 0:
            raise ValueError("ContourData length is not divisible by 3.")
        points = data.reshape(-1, 3)

        expected_points = getattr(contour, "NumberOfContourPoints", None)
        if expected_points is not None and int(expected_points) != points.shape[0]:
            raise ValueError(
                "NumberOfContourPoints does not match ContourData length."
            )

        relative = points - image_position
        plane_offsets = relative @ slice_normal
        if np.max(np.abs(plane_offsets)) > contour_plane_tolerance_mm:
            raise ValueError(
                f"Contour points for ROI {roi_name!r} are not coplanar with "
                f"their referenced CT slice within {contour_plane_tolerance_mm} mm."
            )

        column_coordinates = (
            relative @ column_index_direction
        ) / column_spacing
        row_coordinates = (
            relative @ row_index_direction
        ) / row_spacing

        polygon = _polygon_mask(
            ct.geometry.rows,
            ct.geometry.columns,
            row_coordinates,
            column_coordinates,
        )

        if geometric_type == "CLOSEDPLANAR_XOR":
            mask[slice_index] ^= polygon
        else:
            mask[slice_index] |= polygon

        contour_count += 1
        referenced_slices.add(slice_index)

    if len(geometric_types) > 1:
        raise ValueError(
            f"ROI {roi_name!r} mixes CLOSED_PLANAR and CLOSEDPLANAR_XOR contours."
        )

    if not np.any(mask):
        raise ValueError(f"Rasterized ROI {roi_name!r} contains no CT voxels.")

    return RoiMask(
        roi_number=roi_number,
        roi_name=roi_name,
        mask=mask,
        contour_count=contour_count,
        referenced_slice_count=len(referenced_slices),
    )
