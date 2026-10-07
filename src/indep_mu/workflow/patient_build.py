from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path

import numpy as np

from indep_mu.dicom.case import discover_dicom_case
from indep_mu.dicom.ct import load_ct_series
from indep_mu.dicom.rtplan import RtPlan, load_rtplan
from indep_mu.dicom.rtstruct import build_roi_mask
from indep_mu.patient_model.discretize import discretize_material_blends
from indep_mu.patient_model.egsphant import write_egsphant
from indep_mu.patient_model.external_mask import derive_external_mask_from_ct
from indep_mu.patient_model.hu_red import DRT120KV
from indep_mu.patient_model.model import build_patient_model
from indep_mu.patient_model.pegsless import (
    build_pegsless_media_set,
    write_media_definition,
)


@dataclass(frozen=True)
class PatientArtifactBuild:
    egsphant_path: Path
    media_definition_path: Path
    summary_path: Path
    medium_count: int
    patient_volume_cm3: float
    low_hu_count: int
    high_hu_count: int


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _prepare_outputs(
    output_directory: Path,
    *,
    overwrite: bool,
) -> tuple[Path, Path, Path]:
    output_directory.mkdir(parents=True, exist_ok=True)
    egsphant = output_directory / "patient.egsphant"
    media = output_directory / "media.pegsless"
    summary = output_directory / "patient_model_summary.json"

    existing = [path for path in (egsphant, media, summary) if path.exists()]
    if existing and not overwrite:
        names = ", ".join(path.name for path in existing)
        raise FileExistsError(
            f"Patient-model output already exists: {names}. "
            "Use overwrite=True only after reviewing the previous result."
        )
    return egsphant, media, summary


def _common_plan_isocenter(
    plan: RtPlan,
    *,
    tolerance_mm: float = 1e-3,
) -> tuple[float, float, float]:
    values: list[np.ndarray] = []

    for beam in plan.beams:
        if not beam.control_points:
            raise ValueError(f"Beam {beam.number} has no control points.")
        value = beam.control_points[0].isocenter_position_mm
        if value is None:
            raise ValueError(
                f"Beam {beam.number} has no inherited IsocenterPosition."
            )
        values.append(np.asarray(value, dtype=np.float64))

    reference = values[0]
    for beam, value in zip(plan.beams[1:], values[1:]):
        if not np.allclose(
            value,
            reference,
            rtol=0.0,
            atol=tolerance_mm,
        ):
            raise ValueError(
                "CT-derived external mask requires a common plan isocenter; "
                f"beam {beam.number} differs by more than {tolerance_mm} mm."
            )

    return tuple(float(value) for value in reference)


def build_patient_artifacts(
    dicom_directory: str | Path,
    output_directory: str | Path,
    *,
    patient_roi_name: str | None = None,
    derive_external_from_isocenter: bool = False,
    external_threshold_hu: float = -500.0,
    external_closing_iterations: int = 1,
    mixture_bins: int,
    out_of_range: str = "raise",
    overwrite: bool = False,
) -> PatientArtifactBuild:
    """Build a de-identified EGSnrc patient phantom from one DICOM case.

    The CT series is selected through the RTSTRUCT reference chain, never by
    directory order or instance count.

    Patient-mask source is always explicit:
    - a named RTSTRUCT external/patient ROI; or
    - research-only CT segmentation seeded by the common treatment isocenter.

    The CT-derived path is never selected implicitly.  It is intended for
    validation cases where Monaco did not export an external contour.
    """

    has_roi = patient_roi_name is not None and bool(patient_roi_name.strip())
    if has_roi == bool(derive_external_from_isocenter):
        raise ValueError(
            "Choose exactly one patient-mask source: patient_roi_name or "
            "derive_external_from_isocenter=True."
        )
    if mixture_bins < 1:
        raise ValueError("mixture_bins must be >= 1.")
    if out_of_range not in {"raise", "clip"}:
        raise ValueError("out_of_range must be 'raise' or 'clip'.")

    case = discover_dicom_case(dicom_directory)
    if not case.passed:
        details = "; ".join(f"{m.code}: {m.message}" for m in case.errors)
        raise ValueError(f"DICOM case preflight failed: {details}")

    ct = load_ct_series(
        dicom_directory,
        series_instance_uid=case.ct_series.series_instance_uid,
    )

    patient_roi_summary: dict[str, object] | None = None
    external_mask_summary: dict[str, object] | None = None

    if has_roi:
        assert patient_roi_name is not None
        roi = build_roi_mask(
            ct,
            case.rtstruct.path,
            patient_roi_name,
        )
        patient_mask = roi.mask
        patient_roi_summary = {
            "name": roi.roi_name,
            "number": roi.roi_number,
            "contour_count": roi.contour_count,
            "referenced_slice_count": roi.referenced_slice_count,
        }
        mask_source = "RTSTRUCT_ROI"
    else:
        plan = load_rtplan(case.rtplan.path)
        isocenter = _common_plan_isocenter(plan)
        derived = derive_external_mask_from_ct(
            ct,
            seed_patient_mm=isocenter,
            threshold_hu=external_threshold_hu,
            closing_iterations=external_closing_iterations,
            fill_holes_per_slice=True,
        )
        patient_mask = derived.mask
        mask_source = "CT_DERIVED_RESEARCH"
        diagnostics = derived.diagnostics
        external_mask_summary = {
            "seed_patient_mm": list(diagnostics.seed_patient_mm),
            "threshold_hu": diagnostics.threshold_hu,
            "closing_iterations": int(external_closing_iterations),
            "connected_component_count": diagnostics.connected_component_count,
            "selected_component_label": diagnostics.selected_component_label,
            "selected_component_voxel_count_before_fill": (
                diagnostics.selected_component_voxel_count_before_fill
            ),
            "final_voxel_count": diagnostics.final_voxel_count,
            "touches_ct_border": diagnostics.touches_ct_border,
            "validation_status": "RESEARCH_ONLY",
        }

    model = build_patient_model(
        ct,
        patient_mask,
        calibration=DRT120KV,
        out_of_range=out_of_range,
    )
    discrete = discretize_material_blends(
        model.material_blends,
        mixture_bins=mixture_bins,
    )
    media_set = build_pegsless_media_set(discrete)

    output_root = Path(output_directory)
    egsphant_path, media_path, summary_path = _prepare_outputs(
        output_root,
        overwrite=overwrite,
    )

    geometry = write_egsphant(
        egsphant_path,
        ct=ct,
        model=model,
        discrete=discrete,
    )
    write_media_definition(media_path, media_set)

    body_hu = np.asarray(model.hu[model.patient_mask], dtype=np.float64)
    body_density = np.asarray(
        model.mass_density_g_cm3[model.patient_mask],
        dtype=np.float64,
    )

    if patient_roi_summary is not None:
        patient_roi_summary = {
            **patient_roi_summary,
            "voxel_count": model.patient_voxel_count,
            "volume_cm3": model.patient_volume_cm3,
        }

    summary = {
        "schema_version": 2,
        "source": {
            "ct_series_instance_uid": ct.geometry.series_instance_uid,
            "frame_of_reference_uid": ct.geometry.frame_of_reference_uid,
            "rtstruct_sop_instance_uid": case.rtstruct.sop_instance_uid,
            "rtplan_sop_instance_uid": case.rtplan.sop_instance_uid,
        },
        "ct": {
            "shape_zyx": [int(value) for value in ct.hu.shape],
            "pixel_spacing_mm": [
                float(value) for value in ct.geometry.pixel_spacing_mm
            ],
            "slice_spacing_mm": float(ct.geometry.slice_spacing_mm),
            "patient_position": ct.geometry.patient_position,
            "kvp": ct.geometry.kvp,
            "calibration": "DICOM3.DRT120kV",
            "calibration_hu_range": [DRT120KV.hu_min, DRT120KV.hu_max],
        },
        "patient_mask": {
            "source_type": mask_source,
            "voxel_count": model.patient_voxel_count,
            "volume_cm3": model.patient_volume_cm3,
            "rtstruct_roi": patient_roi_summary,
            "ct_derived": external_mask_summary,
        },
        # Kept for compatibility with early validation scripts.
        "patient_roi": patient_roi_summary,
        "patient_model": {
            "out_of_range_policy": out_of_range,
            "low_hu_count": model.low_hu_count,
            "high_hu_count": model.high_hu_count,
            "body_hu_min": float(np.min(body_hu)),
            "body_hu_max": float(np.max(body_hu)),
            "body_density_min_g_cm3": float(np.min(body_density)),
            "body_density_max_g_cm3": float(np.max(body_density)),
            "cropped_shape_zyx": [int(value) for value in model.shape],
            "mixture_bins": int(mixture_bins),
            "medium_count": len(discrete.media),
            "media": [
                {
                    "name": item.egsphant_name,
                    "material_a": item.material_a,
                    "material_b": item.material_b,
                    "fraction_b": item.fraction_b,
                }
                for item in discrete.media
            ],
        },
        "egsphant_geometry": {
            "shape_zyx": [int(value) for value in geometry.shape_zyx],
            "x_direction_patient": list(geometry.x_direction_patient),
            "y_direction_patient": list(geometry.y_direction_patient),
            "z_direction_patient": list(geometry.z_direction_patient),
            "x_bound_cm": [
                float(geometry.x_bound_cm[0]),
                float(geometry.x_bound_cm[-1]),
            ],
            "y_bound_cm": [
                float(geometry.y_bound_cm[0]),
                float(geometry.y_bound_cm[-1]),
            ],
            "z_bound_cm": [
                float(geometry.z_bound_cm[0]),
                float(geometry.z_bound_cm[-1]),
            ],
        },
        "artifacts": {
            "egsphant": {
                "name": egsphant_path.name,
                "sha256": _sha256(egsphant_path),
            },
            "media_definition": {
                "name": media_path.name,
                "sha256": _sha256(media_path),
            },
        },
        "privacy": {
            "patient_name_stored": False,
            "patient_id_stored": False,
        },
    }

    summary_path.write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    return PatientArtifactBuild(
        egsphant_path=egsphant_path,
        media_definition_path=media_path,
        summary_path=summary_path,
        medium_count=len(discrete.media),
        patient_volume_cm3=model.patient_volume_cm3,
        low_hu_count=model.low_hu_count,
        high_hu_count=model.high_hu_count,
    )
