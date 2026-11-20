"""Research-only patient-mask quality checks; not a couch classifier.

For axial HFS, increasing DICOM row index follows +patient Y (posterior).
The posterior survey measures thresholded material behind the selected
external-mask silhouette. It is a review aid, NOT proof that support has been
excluded: a connected couch may already be inside the selected mask.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import ndimage

from indep_mu.dicom.ct import CtSeries


@dataclass(frozen=True)
class MaskQaReport:
    voxel_count: int
    volume_cm3: float
    touches_ct_border: bool
    calibration_low_count: int
    calibration_high_count: int
    calibration_low_fraction: float
    calibration_high_fraction: float
    posterior_survey_voxel_count: int
    posterior_survey_slice_count: int
    posterior_survey_slice_indices: tuple[int, ...]
    posterior_island_voxel_count: int = 0
    posterior_island_slice_indices: tuple[int, ...] = ()
    support_review_status: str = "MANUAL_REVIEW_REQUIRED"

    @property
    def safe_to_assert_support_absent(self) -> bool:
        # A CT threshold/silhouette analysis cannot verify absence of a
        # connected couch or immobilization board without an independent
        # body/support segmentation reference.
        return False


@dataclass(frozen=True)
class MaskComparison:
    baseline_voxel_count: int
    candidate_voxel_count: int
    dice: float
    added_voxel_count: int
    removed_voxel_count: int
    relative_volume_change_percent: float


def _is_axial_hfs(ct: CtSeries) -> bool:
    expected = np.asarray([1, 0, 0, 0, 1, 0], dtype=np.float64)
    observed = np.asarray(ct.geometry.image_orientation_patient, dtype=np.float64)
    return (
        (ct.geometry.patient_position or "").upper() == "HFS"
        and np.allclose(observed, expected, rtol=0.0, atol=1e-5)
    )


@dataclass(frozen=True)
class PosteriorIslandSurvey:
    suspicious_voxel_count: int
    suspicious_slice_indices: tuple[int, ...]


def survey_posterior_mask_islands(
    mask: np.ndarray,
    *,
    min_component_voxels: int = 20,
    min_component_width_voxels: int = 15,
    gap_rows: int = 2,
) -> PosteriorIslandSurvey:
    """Flag detached posterior components already INCLUDED in a 3-D mask.

    A mask selected as one 3-D connected component can be split into several
    2-D components on many axial slices. This survey detects posterior
    horizontal islands that are *inside* that 3-D component; the previous
    outside-mask posterior survey could not detect these.

    The largest per-slice 2-D component is a geometric reference only.
    Findings are candidates, not confirmed couch pixels. A zero result
    cannot clear an attached support object.
    """

    selected = np.asarray(mask, dtype=bool)
    if selected.ndim != 3 or not selected.any():
        raise ValueError("A non-empty 3-D mask is required.")
    if min_component_voxels < 1 or min_component_width_voxels < 1:
        raise ValueError("Component thresholds must be >= 1.")
    if gap_rows < 0:
        raise ValueError("gap_rows must be non-negative.")

    structure = ndimage.generate_binary_structure(2, 2)
    suspicious_voxels = 0
    suspicious_slices: list[int] = []
    for z_index in range(selected.shape[0]):
        plane = selected[z_index]
        labels, number = ndimage.label(plane, structure=structure)
        if number < 2:
            continue
        sizes = np.bincount(labels.ravel())
        sizes[0] = 0
        reference_label = int(np.argmax(sizes))
        objects = ndimage.find_objects(labels)
        reference = objects[reference_label - 1]
        if reference is None:
            continue
        reference_posterior_stop = int(reference[0].stop)

        found = 0
        for label, slices in enumerate(objects, start=1):
            if slices is None or label == reference_label:
                continue
            if int(sizes[label]) < min_component_voxels:
                continue
            if int(slices[1].stop - slices[1].start) < min_component_width_voxels:
                continue
            if int(slices[0].start) < reference_posterior_stop + gap_rows:
                continue
            found += int(sizes[label])

        if found:
            suspicious_voxels += found
            suspicious_slices.append(z_index)

    return PosteriorIslandSurvey(
        suspicious_voxel_count=suspicious_voxels,
        suspicious_slice_indices=tuple(suspicious_slices),
    )


def review_patient_mask(
    ct: CtSeries,
    mask: np.ndarray,
    *,
    hu_min: float = -1000.0,
    hu_max: float = 2009.0,
    posterior_threshold_hu: float = -500.0,
    posterior_search_mm: float = 50.0,
    posterior_gap_mm: float = 2.0,
) -> MaskQaReport:
    """Summarize out-of-range HU and posterior material for axial HFS CT.

    A positive posterior survey flags material behind the mask that merits
    visual inspection. A zero result does *not* validate a connected couch.
    """

    if not _is_axial_hfs(ct):
        raise ValueError("Posterior survey requires axial HFS CT geometry.")
    if not (np.isfinite(hu_min) and np.isfinite(hu_max) and hu_min < hu_max):
        raise ValueError("Invalid HU calibration range.")
    if not np.isfinite(posterior_threshold_hu):
        raise ValueError("Posterior HU threshold must be finite.")
    if not (0 <= posterior_gap_mm < posterior_search_mm):
        raise ValueError("Posterior survey requires 0 <= gap < search distance.")

    selected = np.asarray(mask, dtype=bool)
    hu = np.asarray(ct.hu)
    if selected.shape != hu.shape or hu.ndim != 3:
        raise ValueError("Mask must match the CT (z, y, x) volume.")
    count = int(np.count_nonzero(selected))
    if count == 0:
        raise ValueError("Patient mask is empty.")

    low = int(np.count_nonzero(selected & (hu < hu_min)))
    high = int(np.count_nonzero(selected & (hu > hu_max)))
    touches = bool(
        selected[0].any()
        or selected[-1].any()
        or selected[:, 0, :].any()
        or selected[:, -1, :].any()
        or selected[:, :, 0].any()
        or selected[:, :, -1].any()
    )

    # The survey is restricted to the lateral span occupied by the body
    # on each slice. Search begins behind the posterior-most mask voxel.
    # This catches separate table/support objects but cannot detect a
    # support component already merged into the body mask.
    row_spacing_mm = float(ct.geometry.pixel_spacing_mm[0])
    first_row_offset = max(1, int(np.ceil(posterior_gap_mm / row_spacing_mm)))
    last_row_offset = int(np.ceil(posterior_search_mm / row_spacing_mm))
    total_posterior_voxels = 0
    suspicious_slices: list[int] = []

    for z_index in range(selected.shape[0]):
        plane = selected[z_index]
        if not plane.any():
            continue
        rows = np.flatnonzero(np.any(plane, axis=1))
        cols = np.flatnonzero(np.any(plane, axis=0))
        last_body_row = int(rows[-1])
        start = last_body_row + first_row_offset
        stop = min(last_body_row + last_row_offset + 1, plane.shape[0])
        if start >= stop:
            continue
        xstart, xstop = int(cols[0]), int(cols[-1]) + 1
        behind = hu[z_index, start:stop, xstart:xstop]
        outside = ~selected[z_index, start:stop, xstart:xstop]
        found = int(
            np.count_nonzero(
                outside & (behind > posterior_threshold_hu)
            )
        )
        total_posterior_voxels += found
        if found:
            suspicious_slices.append(z_index)

    posterior_islands = survey_posterior_mask_islands(selected)

    return MaskQaReport(
        voxel_count=count,
        volume_cm3=count * ct.geometry.voxel_volume_mm3 / 1000.0,
        touches_ct_border=touches,
        calibration_low_count=low,
        calibration_high_count=high,
        calibration_low_fraction=low / count,
        calibration_high_fraction=high / count,
        posterior_survey_voxel_count=total_posterior_voxels,
        posterior_survey_slice_count=len(suspicious_slices),
        posterior_survey_slice_indices=tuple(suspicious_slices),
        posterior_island_voxel_count=posterior_islands.suspicious_voxel_count,
        posterior_island_slice_indices=posterior_islands.suspicious_slice_indices,
    )


def compare_external_masks(
    baseline: np.ndarray,
    candidate: np.ndarray,
) -> MaskComparison:
    """Compute exact Dice and additions/removals without dose assumptions."""

    a = np.asarray(baseline, dtype=bool)
    b = np.asarray(candidate, dtype=bool)
    if a.shape != b.shape or a.ndim != 3:
        raise ValueError("Both masks must have the same 3D shape.")

    count_a = int(np.count_nonzero(a))
    count_b = int(np.count_nonzero(b))
    if count_a == 0 or count_b == 0:
        raise ValueError("Both masks must be non-empty.")

    shared = int(np.count_nonzero(a & b))
    added = int(np.count_nonzero(b & ~a))
    removed = int(np.count_nonzero(a & ~b))
    return MaskComparison(
        baseline_voxel_count=count_a,
        candidate_voxel_count=count_b,
        dice=2.0 * shared / (count_a + count_b),
        added_voxel_count=added,
        removed_voxel_count=removed,
        relative_volume_change_percent=100.0 * (count_b - count_a) / count_a,
    )
