from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from indep_mu.dicom.case import discover_dicom_case
from indep_mu.dicom.ct import load_ct_series
from indep_mu.dicom.rtplan import RtPlan, load_rtplan
from indep_mu.patient_model.external_mask import derive_external_mask_from_ct
from indep_mu.patient_model.hu_red import DRT120KV
from indep_mu.patient_model.model import build_patient_model


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

    reference = points[0]
    for beam, point in zip(plan.beams[1:], points[1:]):
        if not np.allclose(point, reference, rtol=0.0, atol=tolerance_mm):
            raise ValueError(
                "CT-derived patient diagnostics require a common isocenter; "
                f"beam {beam.number} differs by more than {tolerance_mm} mm."
            )

    return tuple(float(value) for value in reference)


def _sample_ct_hu_at_patient_point(
    ct,
    point_patient_mm: tuple[float, float, float],
) -> float:
    from indep_mu.patient_model.external_mask import (
        patient_point_to_ct_index_zyx,
    )

    index = patient_point_to_ct_index_zyx(ct, point_patient_mm)
    return float(ct.hu[index])


def analyze_ct_derived_patient(
    dicom_directory: str | Path,
    *,
    threshold_hu: float = -500.0,
    closing_iterations: int = 1,
    max_seed_snap_distance_mm: float = 50.0,
) -> CtDerivedPatientDiagnostics:
    """Analyze a research CT-derived external mask without creating .egsphant.

    HU outside the measured scanner calibration are clipped only for the
    purpose of deriving density statistics and are always counted explicitly.
    """

    case = discover_dicom_case(dicom_directory)
    if not case.passed:
        details = "; ".join(f"{m.code}: {m.message}" for m in case.errors)
        raise ValueError(f"DICOM case preflight failed: {details}")

    ct = load_ct_series(
        dicom_directory,
        series_instance_uid=case.ct_series.series_instance_uid,
    )
    plan = load_rtplan(case.rtplan.path)
    isocenter = _common_plan_isocenter(plan)

    derived = derive_external_mask_from_ct(
        ct,
        seed_patient_mm=isocenter,
        threshold_hu=threshold_hu,
        closing_iterations=closing_iterations,
        fill_holes_per_slice=True,
        max_seed_snap_distance_mm=max_seed_snap_distance_mm,
    )

    model = build_patient_model(
        ct,
        derived.mask,
        calibration=DRT120KV,
        out_of_range="clip",
    )

    body_hu = np.asarray(model.hu[model.patient_mask], dtype=np.float64)
    body_density = np.asarray(
        model.mass_density_g_cm3[model.patient_mask],
        dtype=np.float64,
    )
    diagnostics = derived.diagnostics

    return CtDerivedPatientDiagnostics(
        threshold_hu=float(threshold_hu),
        closing_iterations=int(closing_iterations),
        seed_patient_mm=isocenter,
        seed_hu=_sample_ct_hu_at_patient_point(ct, isocenter),
        seed_was_in_threshold_component=(
            diagnostics.seed_was_in_threshold_component
        ),
        seed_to_selected_component_distance_mm=(
            diagnostics.seed_to_selected_component_distance_mm
        ),
        connected_component_count=diagnostics.connected_component_count,
        selected_component_voxel_count_before_fill=(
            diagnostics.selected_component_voxel_count_before_fill
        ),
        final_voxel_count=diagnostics.final_voxel_count,
        touches_ct_border=diagnostics.touches_ct_border,
        volume_cm3=model.patient_volume_cm3,
        body_hu_min=float(np.min(body_hu)),
        body_hu_max=float(np.max(body_hu)),
        below_calibration_count=model.low_hu_count,
        above_calibration_count=model.high_hu_count,
        density_min_g_cm3=float(np.min(body_density)),
        density_max_g_cm3=float(np.max(body_density)),
    )


def analyze_external_threshold_sensitivity(
    dicom_directory: str | Path,
    *,
    thresholds_hu: tuple[float, ...] = (-600.0, -500.0, -400.0),
    closing_iterations: int = 1,
) -> tuple[CtDerivedPatientDiagnostics, ...]:
    if not thresholds_hu:
        raise ValueError("At least one external-mask threshold is required.")
    return tuple(
        analyze_ct_derived_patient(
            dicom_directory,
            threshold_hu=value,
            closing_iterations=closing_iterations,
        )
        for value in thresholds_hu
    )
