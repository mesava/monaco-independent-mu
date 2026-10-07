from __future__ import annotations

import json
from pathlib import Path

import click

from indep_mu.dicom.case import discover_dicom_case
from indep_mu.dicom.ct import load_ct_series
from indep_mu.dicom.plan_validation import run_plan_preflight
from indep_mu.dicom.rtplan import load_rtplan
from indep_mu.dicom.rtstruct import list_rtstruct_rois
from indep_mu.montecarlo.egsnrc_backend import EgsnrcReferenceBackend
from indep_mu.patient_model.hu_red import DRT120KV
from indep_mu.workflow.fingerprint import build_technical_case_fingerprint
from indep_mu.workflow.patient_build import build_patient_artifacts
from indep_mu.workflow.preflight import run_transport_preflight


@click.group()
def main() -> None:
    """Независимая проверка MU и дозы для Monaco."""


@main.command("mc-env")
@click.option("--beamnrc", default="beamnrc", show_default=True)
@click.option("--dosxyznrc", default="dosxyznrc", show_default=True)
@click.option(
    "--hen-house",
    type=click.Path(file_okay=False, path_type=Path),
)
@click.option(
    "--egs-home",
    type=click.Path(file_okay=False, path_type=Path),
)
def check_mc_environment(
    beamnrc: str,
    dosxyznrc: str,
    hen_house: Path | None,
    egs_home: Path | None,
) -> None:
    """Проверить доступность reference EGSnrc backend."""

    backend = EgsnrcReferenceBackend(
        beamnrc_executable=beamnrc,
        dosxyznrc_executable=dosxyznrc,
        hen_house=hen_house,
        egs_home=egs_home,
    )
    errors = backend.validate_environment()

    click.echo(f"backend:        {backend.capabilities.backend_id}")
    click.echo(
        "BEAMnrc:        "
        + (
            str(backend.resolved_beamnrc())
            if backend.resolved_beamnrc() is not None
            else "<not found>"
        )
    )
    click.echo(
        "DOSXYZnrc:      "
        + (
            str(backend.resolved_dosxyznrc())
            if backend.resolved_dosxyznrc() is not None
            else "<not found>"
        )
    )
    click.echo(f"GPU:            {backend.capabilities.supports_gpu}")

    if errors:
        click.echo("")
        for item in errors:
            click.echo(f"ERROR: {item}")
        raise click.ClickException("EGSnrc environment is not ready.")

    click.echo("status:         READY")


@main.command("rois")
@click.argument(
    "dicom_directory",
    type=click.Path(exists=True, file_okay=False, path_type=Path),
)
def list_rois(dicom_directory: Path) -> None:
    """Показать ROI из RTSTRUCT без вывода идентификаторов пациента."""

    try:
        case = discover_dicom_case(dicom_directory)
        rois = list_rtstruct_rois(case.rtstruct.path)
    except (ValueError, OSError) as exc:
        raise click.ClickException(str(exc)) from exc

    if not rois:
        click.echo("RTSTRUCT does not contain ROI definitions.")
        return

    for roi in rois:
        click.echo(f"{roi.roi_number:4d}  {roi.roi_name}")


@main.command("fingerprint")
@click.argument(
    "dicom_directory",
    type=click.Path(exists=True, file_okay=False, path_type=Path),
)
@click.option(
    "--output",
    type=click.Path(dir_okay=False, path_type=Path),
)
@click.option(
    "--include-roi-names",
    is_flag=True,
    help="Include RTSTRUCT ROI names in the otherwise de-identified fingerprint.",
)
def fingerprint_case(
    dicom_directory: Path,
    output: Path | None,
    include_roi_names: bool,
) -> None:
    """Сформировать обезличенный технический fingerprint DICOM case."""

    try:
        payload = build_technical_case_fingerprint(
            dicom_directory,
            include_roi_names=include_roi_names,
        )
    except (ValueError, OSError) as exc:
        raise click.ClickException(str(exc)) from exc

    serialized = json.dumps(payload, indent=2, ensure_ascii=False) + "\n"
    if output is None:
        click.echo(serialized, nl=False)
        return

    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(serialized, encoding="utf-8")
    click.echo(str(output))


@main.command("build-patient")
@click.argument(
    "dicom_directory",
    type=click.Path(exists=True, file_okay=False, path_type=Path),
)
@click.argument(
    "output_directory",
    type=click.Path(file_okay=False, path_type=Path),
)
@click.option(
    "--patient-roi",
    required=False,
    help="Exact RTSTRUCT ROI name used as the external patient contour.",
)
@click.option(
    "--derive-external-from-isocenter",
    is_flag=True,
    help="Research-only CT-derived external mask seeded by common plan isocenter.",
)
@click.option(
    "--external-threshold-hu",
    type=float,
    default=-500.0,
    show_default=True,
)
@click.option(
    "--external-closing-iterations",
    type=click.IntRange(min=0),
    default=1,
    show_default=True,
)
@click.option(
    "--mixture-bins",
    required=True,
    type=click.IntRange(min=1),
    help="Material-mixture quantisation; must come from sensitivity validation.",
)
@click.option(
    "--out-of-range",
    type=click.Choice(["raise", "clip"], case_sensitive=False),
    default="raise",
    show_default=True,
    help="Policy for patient voxels outside the measured HU calibration.",
)
@click.option(
    "--overwrite",
    is_flag=True,
    help="Replace an existing patient-model artifact set.",
)
def build_patient(
    dicom_directory: Path,
    output_directory: Path,
    patient_roi: str | None,
    derive_external_from_isocenter: bool,
    external_threshold_hu: float,
    external_closing_iterations: int,
    mixture_bins: int,
    out_of_range: str,
    overwrite: bool,
) -> None:
    """Построить .egsphant и PEGSless media из DICOM case.

    Источник external mask всегда выбирается явно: RTSTRUCT ROI либо
    research-only CT-derived segmentation.
    """

    try:
        result = build_patient_artifacts(
            dicom_directory,
            output_directory,
            patient_roi_name=patient_roi,
            derive_external_from_isocenter=derive_external_from_isocenter,
            external_threshold_hu=external_threshold_hu,
            external_closing_iterations=external_closing_iterations,
            mixture_bins=mixture_bins,
            out_of_range=out_of_range.lower(),
            overwrite=overwrite,
        )
    except (ValueError, OSError) as exc:
        raise click.ClickException(str(exc)) from exc

    click.echo(f"egsphant:      {result.egsphant_path}")
    click.echo(f"media:         {result.media_definition_path}")
    click.echo(f"summary:       {result.summary_path}")
    click.echo(f"media count:   {result.medium_count}")
    click.echo(f"patient volume:{result.patient_volume_cm3:.3f} cm3")
    click.echo(
        "HU outside LUT: "
        f"low={result.low_hu_count}, high={result.high_hu_count}"
    )


@main.command("inspect")
@click.argument(
    "dicom_directory",
    type=click.Path(exists=True, file_okay=False, path_type=Path),
)
def inspect_case(dicom_directory: Path) -> None:
    """Проверить DICOM-комплект без выполнения Monte Carlo."""

    try:
        case = discover_dicom_case(dicom_directory)
        ct = load_ct_series(
            dicom_directory,
            series_instance_uid=case.ct_series.series_instance_uid,
        )
        plan = load_rtplan(case.rtplan.path)
        plan_report = run_plan_preflight(plan)
        transport_report = run_transport_preflight(ct, plan)
        calibration_summary = ct.calibration_range_summary(
            DRT120KV.hu_min,
            DRT120KV.hu_max,
        )
    except (ValueError, OSError) as exc:
        raise click.ClickException(str(exc)) from exc

    case_status = "PASS" if case.passed else "FAIL"
    plan_status = "PASS" if plan_report.passed else "FAIL"
    transport_status = "PASS" if transport_report.passed else "FAIL"

    click.echo(f"DICOM case preflight: {case_status}")
    click.echo(f"Plan preflight:       {plan_status}")
    click.echo(f"Transport preflight:  {transport_status}")
    click.echo(f"CT instances:         {len(case.ct_series.instances)}")
    click.echo(
        "CT matrix:            "
        f"{ct.hu.shape[2]} x {ct.hu.shape[1]} x {ct.hu.shape[0]}"
    )
    click.echo(
        "CT spacing [mm]:      "
        f"{ct.geometry.pixel_spacing_mm[1]:g} x "
        f"{ct.geometry.pixel_spacing_mm[0]:g} x "
        f"{ct.geometry.slice_spacing_mm:g}"
    )
    click.echo(
        f"PatientPosition:       {ct.geometry.patient_position or '<missing>'}"
    )
    click.echo(
        "CT KVP:               "
        + (
            f"{ct.geometry.kvp:g} kV"
            if ct.geometry.kvp is not None
            else "<missing>"
        )
    )
    click.echo(
        "HU observed:          "
        f"{calibration_summary['hu_observed_min']:g} .. "
        f"{calibration_summary['hu_observed_max']:g}"
    )
    click.echo(
        "HU outside 120kV LUT: "
        f"low={calibration_summary['below_calibration_count']}, "
        f"high={calibration_summary['above_calibration_count']} "
        "(whole CT, not BODY-restricted)"
    )
    click.echo(f"RTDOSE objects:       {len(case.rtdoses)}")
    click.echo(f"Treatment beams:      {len(plan.beams)}")
    click.echo("")

    for feature in plan_report.beam_features:
        click.echo(
            "Beam "
            f"{feature.beam_number}: "
            f"{feature.delivery_class}, "
            f"{feature.beam_meterset_mu:.3f} MU, "
            f"{feature.control_point_count} CP"
        )

    all_messages = (
        tuple(case.messages)
        + tuple(plan_report.messages)
        + tuple(transport_report.messages)
    )
    if all_messages:
        click.echo("")
        click.echo("Preflight messages:")
        for item in all_messages:
            beam_suffix = (
                f" [beam {item.beam_number}]"
                if hasattr(item, "beam_number") and item.beam_number is not None
                else ""
            )
            cp_suffix = (
                f" [CP {item.control_point_index}]"
                if (
                    hasattr(item, "control_point_index")
                    and item.control_point_index is not None
                )
                else ""
            )
            click.echo(
                f"{item.severity}: {item.code}{beam_suffix}{cp_suffix}: "
                f"{item.message}"
            )

    if not case.passed or not plan_report.passed or not transport_report.passed:
        raise click.ClickException(
            "Preflight failed. Monte Carlo calculation must not start."
        )


if __name__ == "__main__":
    main()
