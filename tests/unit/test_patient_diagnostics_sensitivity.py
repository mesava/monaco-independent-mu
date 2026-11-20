"""Threshold sensitivity is evaluated once per CT, not once per file scan."""

from __future__ import annotations

import numpy as np

from indep_mu.dicom.ct import CtGeometry, CtSeries
import indep_mu.workflow.patient_diagnostics as diagnostics


def _synthetic_ct() -> CtSeries:
    shape = (9, 15, 17)
    hu = np.full(shape, -1000.0, dtype=np.float32)
    hu[2:7, 4:11, 4:13] = 0.0
    hu[4, 6:9, 7:10] = -1000.0  # An internal low-density cavity.
    hu[2:7, 13:14, 5:12] = 250.0  # Separate posterior support.

    geometry = CtGeometry(
        rows=15,
        columns=17,
        pixel_spacing_mm=(2.0, 2.0),
        slice_spacing_mm=2.5,
        image_orientation_patient=(1.0, 0.0, 0.0, 0.0, 1.0, 0.0),
        first_image_position_patient_mm=(0.0, 0.0, 0.0),
        last_image_position_patient_mm=(0.0, 0.0, 20.0),
        patient_position="HFS",
        frame_of_reference_uid="1.2.3",
        series_instance_uid="1.2.3.4",
        kvp=120.0,
    )
    return CtSeries(
        hu=hu,
        geometry=geometry,
        slice_positions_mm=np.arange(9, dtype=float) * 2.5,
        image_positions_patient_mm=np.column_stack(
            (np.zeros(9), np.zeros(9), np.arange(9, dtype=float) * 2.5)
        ),
        sop_instance_uids=tuple(str(i) for i in range(9)),
    )


def test_threshold_sensitivity_reuses_ct_and_reports_support(monkeypatch) -> None:
    ct = _synthetic_ct()
    load_count = 0

    def fake_load(_path):
        nonlocal load_count
        load_count += 1
        return ct, (16.0, 14.0, 10.0)

    monkeypatch.setattr(diagnostics, "_load_case", fake_load)

    rows = diagnostics.analyze_external_threshold_sensitivity(
        "synthetic",
        thresholds_hu=(-600.0, -500.0, -400.0),
        closing_iterations=0,
    )

    assert load_count == 1
    assert [r.threshold_hu for r in rows] == [-600.0, -500.0, -400.0]
    assert all(r.dice_vs_reference == 1.0 for r in rows)
    assert all(r.mask_added_voxels == 0 for r in rows)
    assert all(r.mask_removed_voxels == 0 for r in rows)
    assert all(r.support_review_status == "MANUAL_REVIEW_REQUIRED" for r in rows)
    assert all(r.posterior_survey_voxel_count > 0 for r in rows)


def test_seed_snap_does_not_use_global_distance_transform(monkeypatch) -> None:
    from scipy import ndimage
    from indep_mu.patient_model.external_mask import derive_external_mask_from_ct

    ct = _synthetic_ct()

    def forbidden_global_edt(*args, **kwargs):
        raise AssertionError("Global EDT should never be used for seed snapping")

    monkeypatch.setattr(ndimage, "distance_transform_edt", forbidden_global_edt)

    result = derive_external_mask_from_ct(
        ct,
        seed_patient_mm=(16.0, 14.0, 10.0),
        threshold_hu=-500.0,
        closing_iterations=0,
        max_seed_snap_distance_mm=20.0,
    )

    assert not result.diagnostics.seed_was_in_threshold_component
    assert 0 < result.diagnostics.seed_to_selected_component_distance_mm <= 20
