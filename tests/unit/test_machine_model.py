from pathlib import Path

import pytest

from indep_mu.beam_model.machine import resolve_energy_model
from indep_mu.dicom.rtplan import load_rtplan
from test_rtplan import _write_plan


def test_standard_6mv_maps_to_flattened_model(tmp_path: Path) -> None:
    path = tmp_path / "plan.dcm"
    _write_plan(path)
    beam = load_rtplan(path).beam(1)

    model = resolve_energy_model(beam)

    assert model.model_id == "6MV"
    assert model.filter_mode == "FF"


def test_override_must_match_nominal_energy(tmp_path: Path) -> None:
    path = tmp_path / "plan.dcm"
    _write_plan(path)
    beam = load_rtplan(path).beam(1)

    with pytest.raises(ValueError, match="does not match"):
        resolve_energy_model(beam, override_model_id="10MV")
