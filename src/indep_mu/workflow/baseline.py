from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import yaml

from indep_mu.beam_model.machine import resolve_energy_model
from indep_mu.dicom.case import discover_dicom_case
from indep_mu.dicom.ct import load_ct_series
from indep_mu.dicom.plan_validation import run_plan_preflight
from indep_mu.dicom.rtdose import load_rtdose
from indep_mu.dicom.rtplan import Beam, load_rtplan
from indep_mu.dicom.rtstruct import list_rtstruct_rois


@dataclass(frozen=True)
class BaselineMismatch:
    path: str
    expected: object
    actual: object


@dataclass(frozen=True)
class BaselineValidationReport:
    expected_path: Path
    mismatches: tuple[BaselineMismatch, ...]

    @property
    def passed(self) -> bool:
        return not self.mismatches


def _common_first_cp_isocenter_mm(
    beams: tuple[Beam, ...],
    *,
    tolerance_mm: float = 1e-3,
) -> list[float] | None:
    values: list[np.ndarray] = []
    for beam in beams:
        if not beam.control_points:
            return None
        point = beam.control_points[0].isocenter_position_mm
        if point is None:
            return None
        values.append(np.asarray(point, dtype=np.float64))

    if not values:
        return None

    reference = values[0]
    if any(
        not np.allclose(value, reference, rtol=0.0, atol=tolerance_mm)
        for value in values[1:]
    ):
        return None
    return [float(value) for value in reference]


def _common_device_definition(
    beams: tuple[Beam, ...],
    device_types: set[str],
):
    matches = []
    for beam in beams:
        found = [
            item
            for item in beam.device_definitions
            if item.device_type in device_types
        ]
        if len(found) != 1:
            return None
        matches.append(found[0])

    reference = matches[0]
    for item in matches[1:]:
        if item != reference:
            return None
    return reference


def build_case_baseline_observation(
    dicom_directory: str | Path,
) -> dict[str, Any]:
    """Build a de-identified observation matching config/validation schema."""

    case = discover_dicom_case(dicom_directory)
    ct = load_ct_series(
        dicom_directory,
        series_instance_uid=case.ct_series.series_instance_uid,
    )
    plan = load_rtplan(case.rtplan.path)
    preflight = run_plan_preflight(plan)
    rois = list_rtstruct_rois(case.rtstruct.path)

    feature_by_number = {
        item.beam_number: item
        for item in preflight.beam_features
    }

    energy_ids = {
        resolve_energy_model(beam).model_id
        for beam in plan.beams
    }
    delivery_classes = {
        feature_by_number[beam.number].delivery_class
        for beam in plan.beams
    }

    mlc = _common_device_definition(plan.beams, {"MLCX", "MLCY"})
    jaw = _common_device_definition(plan.beams, {"ASYMY", "Y"})

    if mlc is not None and mlc.leaf_position_boundaries_mm is not None:
        boundaries = np.asarray(
            mlc.leaf_position_boundaries_mm,
            dtype=np.float64,
        )
        widths = np.diff(boundaries)
        mlc_payload: dict[str, Any] | None = {
            "device_type": mlc.device_type,
            "leaf_pairs": mlc.number_of_leaf_jaw_pairs,
            "first_boundary_mm": float(boundaries[0]),
            "last_boundary_mm": float(boundaries[-1]),
            "nominal_leaf_width_mm": float(np.median(widths)),
            "source_to_device_distance_mm": mlc.source_to_device_distance_mm,
        }
    else:
        mlc_payload = None

    jaw_payload = (
        {
            "device_type": jaw.device_type,
            "pair_count": jaw.number_of_leaf_jaw_pairs,
            "source_to_device_distance_mm": jaw.source_to_device_distance_mm,
        }
        if jaw is not None
        else None
    )

    known_external_names = {"external", "body", "patient"}
    roi_names = {item.roi_name.strip().casefold() for item in rois}
    external_present = bool(roi_names & known_external_names)

    beam_rows = []
    for beam in plan.beams:
        cp0 = beam.control_points[0]
        beam_rows.append(
            {
                "beam_number": beam.number,
                "mu": beam.beam_meterset_mu,
                "control_points": len(beam.control_points),
                "gantry_deg": cp0.gantry_angle_deg,
                "collimator_deg": cp0.collimator_angle_deg,
            }
        )

    doses = [load_rtdose(item.path) for item in case.rtdoses]
    per_beam = any(
        dose.dose_summation_type.upper() == "BEAM"
        for dose in doses
    )

    first_dose = doses[0] if len(doses) == 1 else None
    rtdose_payload = {
        "object_count": len(doses),
        "dose_units": first_dose.dose_units if first_dose is not None else None,
        "dose_type": first_dose.dose_type if first_dose is not None else None,
        "dose_summation_type": (
            first_dose.dose_summation_type
            if first_dose is not None
            else None
        ),
        "per_beam_dose_available": per_beam,
        # Standard RTDOSE does not encode dose-to-medium vs dose-to-water.
        "dose_quantity": "UNKNOWN",
    }

    return {
        "schema_version": 1,
        "ct": {
            "slice_count": int(ct.hu.shape[0]),
            "rows": int(ct.hu.shape[1]),
            "columns": int(ct.hu.shape[2]),
            "pixel_spacing_mm": [
                float(value) for value in ct.geometry.pixel_spacing_mm
            ],
            "slice_spacing_mm": float(ct.geometry.slice_spacing_mm),
            "patient_position": ct.geometry.patient_position,
            "kvp": ct.geometry.kvp,
            "image_orientation_patient": [
                float(value)
                for value in ct.geometry.image_orientation_patient
            ],
        },
        "rtstruct": {
            "roi_count": len(rois),
            "external_roi_present": external_present,
        },
        "rtplan": {
            "number_of_fractions_planned": plan.number_of_fractions_planned,
            "beam_count": len(plan.beams),
            "total_mu_per_fraction": float(
                sum(beam.beam_meterset_mu for beam in plan.beams)
            ),
            "expected_energy_model_id": (
                next(iter(energy_ids))
                if len(energy_ids) == 1
                else None
            ),
            "expected_delivery_class": (
                next(iter(delivery_classes))
                if len(delivery_classes) == 1
                else None
            ),
            "common_isocenter_mm": _common_first_cp_isocenter_mm(plan.beams),
            "mlc": mlc_payload,
            "y_jaws": jaw_payload,
            "beams": beam_rows,
        },
        "rtdose": rtdose_payload,
    }


def _compare_subset(
    expected: object,
    actual: object,
    *,
    path: str,
    mismatches: list[BaselineMismatch],
    float_atol: float,
) -> None:
    if isinstance(expected, dict):
        if not isinstance(actual, dict):
            mismatches.append(BaselineMismatch(path, expected, actual))
            return
        for key, expected_value in expected.items():
            if key in {"case_id", "status"}:
                continue
            child_path = f"{path}.{key}" if path else str(key)
            if key not in actual:
                mismatches.append(
                    BaselineMismatch(child_path, expected_value, "<missing>")
                )
                continue
            _compare_subset(
                expected_value,
                actual[key],
                path=child_path,
                mismatches=mismatches,
                float_atol=float_atol,
            )
        return

    if isinstance(expected, list):
        if not isinstance(actual, list) or len(expected) != len(actual):
            mismatches.append(BaselineMismatch(path, expected, actual))
            return
        for index, (expected_value, actual_value) in enumerate(
            zip(expected, actual)
        ):
            _compare_subset(
                expected_value,
                actual_value,
                path=f"{path}[{index}]",
                mismatches=mismatches,
                float_atol=float_atol,
            )
        return

    if (
        isinstance(expected, (int, float))
        and not isinstance(expected, bool)
        and isinstance(actual, (int, float))
        and not isinstance(actual, bool)
    ):
        if not np.isclose(
            float(expected),
            float(actual),
            rtol=0.0,
            atol=float_atol,
        ):
            mismatches.append(BaselineMismatch(path, expected, actual))
        return

    if expected != actual:
        mismatches.append(BaselineMismatch(path, expected, actual))


def validate_case_against_expected(
    dicom_directory: str | Path,
    expected_yaml: str | Path,
    *,
    float_atol: float = 1e-6,
) -> BaselineValidationReport:
    expected_path = Path(expected_yaml)
    raw = yaml.safe_load(expected_path.read_text(encoding="utf-8"))
    if not isinstance(raw, dict):
        raise ValueError("Validation baseline YAML root must be a mapping.")

    observation = build_case_baseline_observation(dicom_directory)

    mismatches: list[BaselineMismatch] = []
    _compare_subset(
        raw,
        observation,
        path="",
        mismatches=mismatches,
        float_atol=float_atol,
    )
    return BaselineValidationReport(
        expected_path=expected_path,
        mismatches=tuple(mismatches),
    )
