from pathlib import Path

import pytest

from indep_mu.dicom.plan_validation import (
    analyze_beam_delivery,
    run_plan_preflight,
)
from indep_mu.dicom.rtplan import load_rtplan
from test_rtplan import _write_plan


def test_delivery_is_detected_as_vmat(tmp_path: Path) -> None:
    path = tmp_path / "plan.dcm"
    _write_plan(path)

    plan = load_rtplan(path)
    features = analyze_beam_delivery(plan.beam(1))

    assert features.delivery_class == "VMAT"
    assert features.mlc_dynamic
    assert features.gantry_dynamic
    assert features.control_point_count == 3
    assert features.segment_count == 2
    assert features.energies_mv == (6.0,)


def test_preflight_passes_supported_synthetic_plan(tmp_path: Path) -> None:
    path = tmp_path / "plan.dcm"
    _write_plan(path)

    result = run_plan_preflight(load_rtplan(path))

    assert result.passed
    assert not result.errors


def test_preflight_rejects_unsupported_modifier(tmp_path: Path) -> None:
    path = tmp_path / "plan.dcm"
    _write_plan(path)

    import pydicom

    ds = pydicom.dcmread(path)
    ds.BeamSequence[0].NumberOfWedges = 1
    ds.save_as(path, enforce_file_format=True)

    result = run_plan_preflight(load_rtplan(path))

    assert not result.passed
    assert any(
        item.code == "UNSUPPORTED_BEAM_MODIFIER"
        for item in result.errors
    )
