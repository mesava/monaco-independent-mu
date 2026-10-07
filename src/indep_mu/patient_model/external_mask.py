from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import ndimage

from indep_mu.dicom.ct import CtSeries


@dataclass(frozen=True)
class CtExternalMaskDiagnostics:
    threshold_hu: float
    connected_component_count: int
    selected_component_label: int
    selected_component_voxel_count_before_fill: int
    final_voxel_count: int
    touches_ct_border: bool
    seed_patient_mm: tuple[float, float, float]
    seed_was_in_threshold_component: bool
    seed_to_selected_component_distance_mm: float


@dataclass(frozen=True)
class CtDerivedExternalMask:
    """Research-only CT-derived external patient mask."""

    mask: np.ndarray
    diagnostics: CtExternalMaskDiagnostics


def patient_point_to_ct_index_zyx(
    ct: CtSeries,
    point_patient_mm: tuple[float, float, float],
) -> tuple[int, int, int]:
    """Map a patient-space point to the nearest CT voxel centre."""

    point = np.asarray(point_patient_mm, dtype=np.float64)
    if point.shape != (3,) or not np.all(np.isfinite(point)):
        raise ValueError("point_patient_mm must contain three finite values.")

    orientation = np.asarray(
        ct.geometry.image_orientation_patient,
        dtype=np.float64,
    )
    column_direction = orientation[:3]
    row_direction = orientation[3:]
    normal = np.cross(column_direction, row_direction)
    normal /= np.linalg.norm(normal)

    row_spacing, column_spacing = ct.geometry.pixel_spacing_mm

    slice_distances = (
        np.asarray(ct.image_positions_patient_mm, dtype=np.float64) - point
    ) @ normal
    z_index = int(np.argmin(np.abs(slice_distances)))

    relative = point - ct.image_positions_patient_mm[z_index]
    column = float(relative @ column_direction) / column_spacing
    row = float(relative @ row_direction) / row_spacing

    row_index = int(np.rint(row))
    column_index = int(np.rint(column))

    if not (
        0 <= z_index < ct.hu.shape[0]
        and 0 <= row_index < ct.hu.shape[1]
        and 0 <= column_index < ct.hu.shape[2]
    ):
        raise ValueError("Seed point lies outside the CT voxel grid.")

    return z_index, row_index, column_index


def _touches_border(mask: np.ndarray) -> bool:
    return bool(
        np.any(mask[0])
        or np.any(mask[-1])
        or np.any(mask[:, 0, :])
        or np.any(mask[:, -1, :])
        or np.any(mask[:, :, 0])
        or np.any(mask[:, :, -1])
    )


def derive_external_mask_from_ct(
    ct: CtSeries,
    *,
    seed_patient_mm: tuple[float, float, float],
    threshold_hu: float = -500.0,
    closing_iterations: int = 1,
    fill_holes_per_slice: bool = True,
    max_seed_snap_distance_mm: float = 50.0,
) -> CtDerivedExternalMask:
    """Derive a research external mask from CT and a plan-space seed.

    The mask is formed from an explicit HU threshold and connected components.
    The component is selected by the supplied patient-space seed, normally the
    common treatment isocenter.

    If the seed is itself in a low-density cavity, the nearest thresholded
    tissue voxel may be used, but only within an explicit maximum distance.
    This permits an isocenter inside lung or another enclosed low-density
    region without silently choosing the largest CT component.

    The result remains research-only. A patient physically connected to the CT
    couch at the chosen threshold can still merge with support hardware and
    therefore requires independent validation against a trusted external
    contour.
    """

    if not np.isfinite(threshold_hu):
        raise ValueError("threshold_hu must be finite.")
    if closing_iterations < 0:
        raise ValueError("closing_iterations must be non-negative.")
    if (
        not np.isfinite(max_seed_snap_distance_mm)
        or max_seed_snap_distance_mm <= 0
    ):
        raise ValueError("max_seed_snap_distance_mm must be finite and positive.")

    hu = np.asarray(ct.hu, dtype=np.float32)
    candidate = hu > float(threshold_hu)

    if closing_iterations:
        structure = ndimage.generate_binary_structure(rank=3, connectivity=1)
        candidate = ndimage.binary_closing(
            candidate,
            structure=structure,
            iterations=closing_iterations,
        )

    structure = ndimage.generate_binary_structure(rank=3, connectivity=1)
    labels, component_count = ndimage.label(candidate, structure=structure)
    if component_count < 1:
        raise ValueError("CT threshold produced no connected components.")

    seed_index = patient_point_to_ct_index_zyx(ct, seed_patient_mm)
    selected_label = int(labels[seed_index])
    seed_was_inside = selected_label != 0
    seed_distance_mm = 0.0

    if selected_label == 0:
        row_spacing, column_spacing = ct.geometry.pixel_spacing_mm
        distances, nearest = ndimage.distance_transform_edt(
            ~candidate,
            sampling=(
                ct.geometry.slice_spacing_mm,
                row_spacing,
                column_spacing,
            ),
            return_indices=True,
        )
        seed_distance_mm = float(distances[seed_index])
        if seed_distance_mm > max_seed_snap_distance_mm:
            raise ValueError(
                "No thresholded patient component is sufficiently close to "
                "the seed point."
            )

        nearest_index = tuple(
            int(nearest[axis][seed_index])
            for axis in range(3)
        )
        selected_label = int(labels[nearest_index])
        if selected_label == 0:
            raise ValueError(
                "Could not resolve a thresholded tissue component near seed."
            )

    selected = labels == selected_label
    before_fill = int(np.count_nonzero(selected))

    if fill_holes_per_slice:
        filled = np.empty_like(selected, dtype=bool)
        for index in range(selected.shape[0]):
            filled[index] = ndimage.binary_fill_holes(selected[index])
        selected = filled

    diagnostics = CtExternalMaskDiagnostics(
        threshold_hu=float(threshold_hu),
        connected_component_count=int(component_count),
        selected_component_label=selected_label,
        selected_component_voxel_count_before_fill=before_fill,
        final_voxel_count=int(np.count_nonzero(selected)),
        touches_ct_border=_touches_border(selected),
        seed_patient_mm=tuple(float(value) for value in seed_patient_mm),
        seed_was_in_threshold_component=seed_was_inside,
        seed_to_selected_component_distance_mm=seed_distance_mm,
    )

    return CtDerivedExternalMask(
        mask=np.asarray(selected, dtype=bool),
        diagnostics=diagnostics,
    )
