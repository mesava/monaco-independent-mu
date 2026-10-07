from __future__ import annotations

from pathlib import Path

import click

from indep_mu.dicom.case import discover_dicom_case
from indep_mu.dicom.plan_validation import run_plan_preflight
from indep_mu.dicom.rtplan import load_rtplan


@click.group()
def main() -> None:
    """Независимая проверка MU и дозы для Monaco."""


@main.command("inspect")
@click.argument(
    "dicom_directory",
    type=click.Path(exists=True, file_okay=False, path_type=Path),
)
def inspect_case(dicom_directory: Path) -> None:
    """Проверить DICOM-комплект без выполнения Monte Carlo."""

    try:
        case = discover_dicom_case(dicom_directory)
        plan = load_rtplan(case.rtplan.path)
        plan_report = run_plan_preflight(plan)
    except (ValueError, OSError) as exc:
        raise click.ClickException(str(exc)) from exc

    case_status = "PASS" if case.passed else "FAIL"
    plan_status = "PASS" if plan_report.passed else "FAIL"

    click.echo(f"DICOM case preflight: {case_status}")
    click.echo(f"Plan preflight:       {plan_status}")
    click.echo(f"CT instances:         {len(case.ct_series.instances)}")
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

    all_messages = tuple(case.messages) + tuple(plan_report.messages)
    if all_messages:
        click.echo("")
        click.echo("Preflight messages:")
        for item in all_messages:
            beam_suffix = (
                f" [beam {item.beam_number}]"
                if hasattr(item, "beam_number") and item.beam_number is not None
                else ""
            )
            click.echo(
                f"{item.severity}: {item.code}{beam_suffix}: {item.message}"
            )

    if not case.passed or not plan_report.passed:
        raise click.ClickException(
            "Preflight failed. Monte Carlo calculation must not start."
        )


if __name__ == "__main__":
    main()
