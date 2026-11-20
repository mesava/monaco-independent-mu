"""Patient-model diagnostics for research CT-derived external masks.

The scanner calibration is evaluated only inside the selected mask. Threshold
sensitivity loads the DICOM case once; this avoids three redundant parses of
the same 187-slice CT and three full patient material-model allocations.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from indep_mu.dicom.case import discover_dicom_case
from indep_mu.dicom.ct import CtSeries, load_ct_series
from indep_mu.dicom.rtplan import RtPlan, load_rtplan
from indep_mu.patient_model.density import red_to_mass_density
from indep_mu.patient_model.external_mask import derive_external_mask_from_ct
from indep_mu.patient_model.hu_red import DRT120KV
from indep_mu.patient_model.mask_qc import (
    MaskQaReport,
    compare_external_masks,
    review_patient_mask,
)


@dataclass(frozen=True)
class CtDerivedPatientDiagnostics:
    threshold_hu: float
    closing_iterations: int
    seed_patient_mm: tuple[float, float, float]
    seed_hu: float
    seed_was_in_threshold_component: bool
    seed_to_selected_component_distance_mm: float
    connected_component_count: int
    selected_component_voxel_count_before_fill: int
    final_voxel_count: int
    touches_ct_border: bool
    volume_cm3: float
    body_hu_min: float
    body_hu_max: float
    below_calibration_count: int
    above_calibration_count: int
    density_min_g_cm3: float
    density_max_g_cm3: float
    posterior_survey_voxel_count: int
    posterior_survey_slice_count: int
    posterior_island_voxel_count: int
    posterior_island_slice_count: int
    support_review_status: str
    reference_threshold_hu: float | None = None
    dice_vs_reference: float | None = None
    mask_added_voxels: int | None = None
    mask_removed_voxels: int | None = None
    relative_volume_change_percent: float | None = None


def _common_plan_isocenter(
    plan: RtPlan,
    *,
    tolerance_mm: float = 1e-3,
) -> tuple[float, float, float]:
    points: list[np.ndarray] = []
    for beam in plan.beams:
        if not beam.control_points:
            raise ValueError(f"Beam {beam.number} has no control points.")
        value = beam.control_points[0].isocenter_position_mm
        if value is None:
            raise ValueError(
                f"Beam {beam.number} has no inherited IsocenterPosition."
            )
        points.append(np.asarray(value, dtype=np.float64))

    if not points:
        raise ValueError("Plan contains no treatment beams.")

    reference = points[0]
    for beam, point in zip(plan.beams[1:], points[1:]):
        if not np.allclose(point, reference, rtol=0.0, atol=tolerance_mm):
            raise ValueError(
                "CT-derived patient diagnostics require a common isocenter; "
                f"beam {beam.number} differs by more than {tolerance_mm} mm."
            )
    return tuple(float(value) for value in reference)


def _load_case(
    dicom_directory: str | Path,
) -> tuple[CtSeries, tuple[float, float, float]]:
    case = discover_dicom_case(dicom_directory)
    if not case.passed:
        details = "; ".join(f"{m.code}: {m.message}" for m in case.errors)
        raise ValueError(f"DICOM case preflight failed: {details}")

    ct = load_ct_series(
        dicom_directory,
        series_instance_uid=case.ct_series.series_instance_uid,
    )
    plan = load_rtplan(case.rtplan.path)
    return ct, _common_plan_isocenter(plan)


def _analyze_loaded(
    ct: CtSeries,
    isocenter: tuple[float, float, float],
    *,
    threshold_hu: float,
    closing_iterations: int,
    max_seed_snap_distance_mm: float,
) -> tuple[CtDerivedPatientDiagnostics, np.ndarray]:
    derived = derive_external_mask_from_ct(
        ct,
        seed_patient_mm=isocenter,
        threshold_hu=threshold_hu,
        closing_iterations=closing_iterations,
        fill_holes_per_slice=True,
        max_seed_snap_distance_mm=max_seed_snap_distance_mm,
    )

    mask: np.ndarray = derived.mask
    qa: MaskQaReport = review_patient_mask(
        ct,
        mask,
        hu_min=DRT120KV.hu_min,
        hu_max=DRT120KV.hu_max,
        posterior_threshold_hu=threshold_hu,
    )

    # Density statistics in the mask only. This matches the model's
    # HU->RED->density conversion but avoids allocating the full material
    # blend field, which is unnecessary for threshold sensitivity.
    body_hu = np.asarray(ct.hu[mask], dtype=np.float64)
    body_red = DRT120KV.to_red(body_hu, out_of_range="clip")
    density = red_to_mass_density(body_red)
    del body_red

    from indep_mu.patient_model.external_mask import patient_point_to_ct_index_zyx

    seed_index = patient_point_to_ct_index_zyx(ct, isocenter)
    source = derived.diagnostics

    return (
        CtDerivedPatientDiagnostics(
            threshold_hu=float(threshold_hu),
            closing_iterations=int(closing_iterations),
            seed_patient_mm=isocenter,
            seed_hu=float(ct.hu[seed_index]),
            seed_was_in_threshold_component=source.seed_was_in_threshold_component,
            seed_to_selected_component_distance_mm=(
                source.seed_to_selected_component_distance_mm
            ),
            connected_component_count=source.connected_component_count,
            selected_component_voxel_count_before_fill=(
                source.selected_component_voxel_count_before_fill
            ),
            final_voxel_count=source.final_voxel_count,
            touches_ct_border=source.touches_ct_border,
            volume_cm3=qa.volume_cm3,
            body_hu_min=float(np.min(body_hu)),
            body_hu_max=float(np.max(body_hu)),
            below_calibration_count=qa.calibration_low_count,
            above_calibration_count=qa.calibration_high_count,
            density_min_g_cm3=float(np.min(density)),
            density_max_g_cm3=float(np.max(density)),
            posterior_survey_voxel_count=qa.posterior_survey_voxel_count,
            posterior_survey_slice_count=qa.posterior_survey_slice_count,
            posterior_island_voxel_count=qa.posterior_island_voxel_count,
            posterior_island_slice_count=len(qa.posterior_island_slice_indices),
            support_review_status=qa.support_review_status,
        ),
        mask,
    )


def analyze_ct_derived_patient(
    dicom_directory: str | Path,
    *,
    threshold_hu: float = -500.0,
    closing_iterations: int = 1,
    max_seed_snap_distance_mm: float = 50.0,
) -> CtDerivedPatientDiagnostics:
    """Analyze one research external mask without creating .egsphant."""

    ct, isocenter = _load_case(dicom_directory)
    result, _ = _analyze_loaded(
        ct,
        isocenter,
        threshold_hu=threshold_hu,
        closing_iterations=closing_iterations,
        max_seed_snap_distance_mm=max_seed_snap_distance_mm,
    )
    return result


def analyze_external_threshold_sensitivity(
    dicom_directory: str | Path,
    *,
    thresholds_hu: tuple[float, ...] = (-600.0, -500.0, -400.0),
    closing_iterations: int = 1,
) -> tuple[CtDerivedPatientDiagnostics, ...]:
    """Compare external masks against -500 HU, or the first chosen threshold.

    All Dice values and added/removed voxel counts are geometric, not dose
    sensitivity measures. Actual support clearance still requires review.
    """

    if not thresholds_hu:
        raise ValueError("At least one external-mask threshold is required.")
    if any(not np.isfinite(x) for x in thresholds_hu):
        raise ValueError("Threshold values must be finite.")

    ct, isocenter = _load_case(dicom_directory)
    reference_threshold = -500.0 if -500.0 in thresholds_hu else thresholds_hu[0]
    baseline, baseline_mask = _analyze_loaded(
        ct,
        isocenter,
        threshold_hu=reference_threshold,
        closing_iterations=closing_iterations,
        max_seed_snap_distance_mm=50.0,
    )

    from dataclasses import replace

    rows = []
    for threshold in thresholds_hu:
        if threshold == reference_threshold:
            result = baseline
            comparison = compare_external_masks(baseline_mask, baseline_mask)
        else:
            result, candidate_mask = _analyze_loaded(
                ct,
                isocenter,
                threshold_hu=threshold,
                closing_iterations=closing_iterations,
                max_seed_snap_distance_mm=50.0,
            )
            comparison = compare_external_masks(baseline_mask, candidate_mask)
            del candidate_mask

        rows.append(
            replace(
                result,
                reference_threshold_hu=float(reference_threshold),
                dice_vs_reference=comparison.dice,
                mask_added_voxels=comparison.added_voxel_count,
                mask_removed_voxels=comparison.removed_voxel_count,
                relative_volume_change_percent=(
                    comparison.relative_volume_change_percent
                ),
            )
        )
    return tuple(rows)
