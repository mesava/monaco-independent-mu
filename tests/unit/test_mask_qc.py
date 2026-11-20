import numpy as np
import pytest

from indep_mu.dicom.ct import CtGeometry, CtSeries
from indep_mu.patient_model.mask_qc import (
    compare_external_masks,
    review_patient_mask,
)


def _ct() -> tuple[CtSeries, np.ndarray]:
    hu = np.full((3, 12, 14), -1000.0, dtype=np.float32)
    mask = np.zeros_like(hu, dtype=bool)
    mask[:, 3:7, 3:10] = True
    hu[mask] = 0.0

    # Distinct posterior couch/support component, not part of the body mask.
    hu[:, 9:11, 4:9] = 200.0
    hu[1, 4, 4] = -1005.0
    hu[1, 5, 5] = 2500.0

    geometry = CtGeometry(
        rows=12,
        columns=14,
        pixel_spacing_mm=(2.0, 2.0),
        slice_spacing_mm=3.0,
        image_orientation_patient=(1.0, 0.0, 0.0, 0.0, 1.0, 0.0),
        first_image_position_patient_mm=(0.0, 0.0, 0.0),
        last_image_position_patient_mm=(0.0, 0.0, 6.0),
        patient_position="HFS",
        frame_of_reference_uid="1.2.3",
        series_instance_uid="1.2.3.4",
        kvp=120.0,
    )
    ct = CtSeries(
        hu=hu,
        geometry=geometry,
        slice_positions_mm=np.asarray([0.0, 3.0, 6.0]),
        image_positions_patient_mm=np.asarray(
            [[0.0, 0.0, 0.0], [0.0, 0.0, 3.0], [0.0, 0.0, 6.0]]
        ),
        sop_instance_uids=("1", "2", "3"),
    )
    return ct, mask


def test_couch_candidate_and_calibration_counts() -> None:
    ct, mask = _ct()
    report = review_patient_mask(
        ct,
        mask,
        posterior_threshold_hu=-500,
        posterior_search_mm=12.0,
        posterior_gap_mm=2.0,
    )

    assert report.voxel_count == 3 * 4 * 7
    assert report.calibration_low_count == 1
    assert report.calibration_high_count == 1
    assert report.calibration_low_fraction == pytest.approx(1 / 84)
    assert report.calibration_high_fraction == pytest.approx(1 / 84)
    assert report.posterior_survey_voxel_count == 3 * 2 * 5
    assert report.posterior_survey_slice_count == 3
    assert report.posterior_survey_slice_indices == (0, 1, 2)
    assert not report.touches_ct_border
    assert report.support_review_status == "MANUAL_REVIEW_REQUIRED"
    assert not report.safe_to_assert_support_absent


def test_absence_of_detached_component_is_not_support_clearance() -> None:
    ct, mask = _ct()
    ct.hu[:, 9:11, 4:9] = -1000

    report = review_patient_mask(ct, mask)
    assert report.posterior_survey_voxel_count == 0
    assert not report.safe_to_assert_support_absent


def test_mask_comparison_counts_and_dice() -> None:
    baseline = np.zeros((2, 3, 4), dtype=bool)
    baseline[0, 1, 1:3] = True
    candidate = baseline.copy()
    candidate[0, 1, 1] = False
    candidate[1, 2, 3] = True

    comparison = compare_external_masks(baseline, candidate)
    assert comparison.baseline_voxel_count == 2
    assert comparison.candidate_voxel_count == 2
    assert comparison.added_voxel_count == 1
    assert comparison.removed_voxel_count == 1
    assert comparison.dice == pytest.approx(0.5)
    assert comparison.relative_volume_change_percent == 0.0


def test_non_hfs_or_oblique_is_not_silently_surveyed() -> None:
    ct, mask = _ct()
    object.__setattr__(ct.geometry, "patient_position", "FFS")

    with pytest.raises(ValueError, match="axial HFS"):
        review_patient_mask(ct, mask)
