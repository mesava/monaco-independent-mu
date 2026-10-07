from dataclasses import replace
from pathlib import Path

from indep_mu.beam_model.agility import validate_agility_dicom_geometry
from indep_mu.dicom.rtplan import (
    BeamLimitingDeviceDefinition,
    load_rtplan,
)
from test_rtplan import _write_plan


def test_synthetic_non_agility_geometry_is_rejected(tmp_path: Path) -> None:
    path = tmp_path / "plan.dcm"
    _write_plan(path)
    beam = load_rtplan(path).beam(1)

    report = validate_agility_dicom_geometry(beam)

    assert not report.passed
    assert any(item.code == "LEAF_PAIR_COUNT" for item in report.errors)


def test_nominal_agility_boundaries_pass_core_checks(tmp_path: Path) -> None:
    path = tmp_path / "plan.dcm"
    _write_plan(path)
    beam = load_rtplan(path).beam(1)

    definitions = tuple(
        BeamLimitingDeviceDefinition(
            device_type=item.device_type,
            number_of_leaf_jaw_pairs=80,
            leaf_position_boundaries_mm=tuple(
                -200.0 + 5.0 * index for index in range(81)
            ),
            source_to_device_distance_mm=349.0,
        )
        if item.device_type == "MLCX"
        else item
        for item in beam.device_definitions
    )
    beam = replace(beam, device_definitions=definitions)

    report = validate_agility_dicom_geometry(beam)

    assert report.passed
    assert report.leaf_pairs == 80
    assert report.nominal_leaf_width_mm == 5.0
    assert report.field_span_mm == 400.0
    assert report.source_to_mlc_distance_mm == 349.0
