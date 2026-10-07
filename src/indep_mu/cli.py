from __future__ import annotations

from pathlib import Path

import click

from indep_mu.dicom.case import discover_dicom_case
from indep_mu.dicom.ct import load_ct_series
from indep_mu.dicom.plan_validation import run_plan_preflight
from indep_mu.dicom.rtplan import load_rtplan
from indep_mu.dicom.rtstruct import list_rtstruct_rois
from indep_mu.patient_model.hu_red import DRT120KV
from indep_mu.workflow.preflight import run_transport_preflight


@click.group()
def main() -> None:
    """Независимая проверка MU и дозы для Monaco."""


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
