from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Any

import numpy as np

from indep_mu.beam_model.agility import validate_agility_dicom_geometry
from indep_mu.beam_model.machine import resolve_energy_model
from indep_mu.dicom.case import discover_dicom_case
from indep_mu.dicom.ct import load_ct_series
from indep_mu.dicom.plan_validation import run_plan_preflight
from indep_mu.dicom.rtdose import load_rtdose
from indep_mu.dicom.rtplan import load_rtplan
from indep_mu.dicom.rtstruct import list_rtstruct_rois
from indep_mu.patient_model.hu_red import DRT120KV
from indep_mu.workflow.preflight import run_transport_preflight


def _hash_identifier(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()[:16]


def build_technical_case_fingerprint(
    dicom_directory: str | Path,
    *,
    include_roi_names: bool = False,
) -> dict[str, Any]:
    """Return a de-identified, reproducible technical DICOM fingerprint.

    PatientName and PatientID are never read into the result. DICOM UIDs are
    represented only by short SHA-256 fingerprints so the output can be stored
    in Git for regression/validation work without copying raw identifiers.
    """

    case = discover_dicom_case(dicom_directory)
    ct = load_ct_series(
        dicom_directory,
        series_instance_uid=case.ct_series.series_instance_uid,
    )
    plan = load_rtplan(case.rtplan.path)
    plan_report = run_plan_preflight(plan)
    transport_report = run_transport_preflight(ct, plan)
    rois = list_rtstruct_rois(case.rtstruct.path)

    calibration = ct.calibration_range_summary(
        DRT120KV.hu_min,
        DRT120KV.hu_max,
    )

    beam_rows: list[dict[str, Any]] = []
    feature_by_number = {
        item.beam_number: item
        for item in plan_report.beam_features
    }

    for beam in plan.beams:
        feature = feature_by_number[beam.number]
        agility = validate_agility_dicom_geometry(beam)
        try:
            energy_model_id = resolve_energy_model(beam).model_id
        except ValueError:
            energy_model_id = None

        cp0 = beam.control_points[0]
        beam_rows.append(
            {
                "beam_number": beam.number,
                "beam_name": beam.name,
                "delivery_class": feature.delivery_class,
                "beam_meterset_mu": beam.beam_meterset_mu,
                "control_point_count": feature.control_point_count,
                "segment_count": feature.segment_count,
                "energy_model_id": energy_model_id,
                "fluence_mode": beam.fluence_mode,
                "fluence_mode_id": beam.fluence_mode_id,
                "gantry_deg_cp0": cp0.gantry_angle_deg,
                "collimator_deg_cp0": cp0.collimator_angle_deg,
                "patient_support_deg_cp0": cp0.patient_support_angle_deg,
                "isocenter_mm_cp0": (
                    list(cp0.isocenter_position_mm)
                    if cp0.isocenter_position_mm is not None
                    else None
                ),
                "agility": {
                    "passed": agility.passed,
                    "mlc_device_type": agility.mlc_device_type,
                    "leaf_pairs": agility.leaf_pairs,
                    "field_span_mm": agility.field_span_mm,
                    "nominal_leaf_width_mm": agility.nominal_leaf_width_mm,
                },
            }
        )

    dose_rows: list[dict[str, Any]] = []
    for dose_ref in case.rtdoses:
        dose = load_rtdose(dose_ref.path)
        row_spacing, column_spacing = dose.geometry.pixel_spacing_mm
        local_offsets = np.asarray(
            dose.geometry.local_frame_offsets_mm,
            dtype=np.float64,
        )
        frame_step_mm = (
            float(np.median(np.abs(np.diff(local_offsets))))
            if local_offsets.size > 1
            else None
        )
        dose_rows.append(
            {
                "sop_uid_hash": _hash_identifier(dose.sop_instance_uid),
                "dose_units": dose.dose_units,
                "dose_type": dose.dose_type,
                "dose_summation_type": dose.dose_summation_type,
                "shape_frc": list(dose.shape),
                "pixel_spacing_mm": [row_spacing, column_spacing],
                "frame_spacing_mm_median": frame_step_mm,
                "dose_min": float(np.min(dose.dose)),
                "dose_max": float(np.max(dose.dose)),
                "referenced_beam_numbers": list(dose.referenced_beam_numbers),
            }
        )

    roi_payload: dict[str, Any] = {
        "count": len(rois),
        "names_included": bool(include_roi_names),
    }
    if include_roi_names:
        roi_payload["names"] = [item.roi_name for item in rois]

    return {
        "schema_version": 1,
        "privacy": {
            "patient_name_included": False,
            "patient_id_included": False,
            "raw_uids_included": False,
        },
        "dicom_case": {
            "passed": case.passed,
            "warnings": [item.code for item in case.warnings],
            "errors": [item.code for item in case.errors],
        },
        "ct": {
            "series_uid_hash": _hash_identifier(
                ct.geometry.series_instance_uid
            ),
            "frame_of_reference_uid_hash": _hash_identifier(
                ct.geometry.frame_of_reference_uid
            ),
            "slice_count": int(ct.hu.shape[0]),
            "rows": int(ct.hu.shape[1]),
            "columns": int(ct.hu.shape[2]),
            "pixel_spacing_mm": list(ct.geometry.pixel_spacing_mm),
            "slice_spacing_mm": ct.geometry.slice_spacing_mm,
            "image_orientation_patient": list(
                ct.geometry.image_orientation_patient
            ),
            "patient_position": ct.geometry.patient_position,
            "kvp": ct.geometry.kvp,
            "hu_observed_min": calibration["hu_observed_min"],
            "hu_observed_max": calibration["hu_observed_max"],
            "below_calibration_count_whole_ct": calibration[
                "below_calibration_count"
            ],
            "above_calibration_count_whole_ct": calibration[
                "above_calibration_count"
            ],
        },
        "rtstruct": roi_payload,
        "rtplan": {
            "plan_label": plan.plan_label,
            "plan_name": plan.plan_name,
            "fraction_group_number": plan.fraction_group_number,
            "number_of_fractions_planned": plan.number_of_fractions_planned,
            "beam_count": len(plan.beams),
            "total_mu_per_fraction": float(
                sum(beam.beam_meterset_mu for beam in plan.beams)
            ),
            "plan_preflight_passed": plan_report.passed,
            "transport_preflight_passed": transport_report.passed,
            "plan_preflight_errors": [
                item.code for item in plan_report.errors
            ],
            "transport_preflight_errors": [
                item.code for item in transport_report.errors
            ],
            "beams": beam_rows,
        },
        "rtdose": dose_rows,
    }
