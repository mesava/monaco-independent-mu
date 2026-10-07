from __future__ import annotations

import os
from pathlib import Path

from indep_mu.montecarlo.egsnrc_backend import EgsnrcReferenceBackend


def _make_executable(path: Path) -> None:
    path.write_text("#!/bin/sh\nexit 0\n", encoding="ascii")
    if os.name != "nt":
        path.chmod(0o755)


def test_egsnrc_backend_capabilities() -> None:
    backend = EgsnrcReferenceBackend()

    assert backend.capabilities.backend_id == "egsnrc-reference"
    assert backend.capabilities.supports_dynamic_mlc
    assert backend.capabilities.supports_dynamic_jaws
    assert backend.capabilities.supports_dynamic_gantry
    assert backend.capabilities.supports_patient_ct
    assert backend.capabilities.supports_dose_to_medium
    assert not backend.capabilities.supports_gpu


def test_environment_validation_with_direct_executables(tmp_path: Path) -> None:
    beamnrc = tmp_path / "beamnrc"
    dosxyz = tmp_path / "dosxyznrc"
    _make_executable(beamnrc)
    _make_executable(dosxyz)

    hen_house = tmp_path / "HEN_HOUSE"
    egs_home = tmp_path / "EGS_HOME"
    hen_house.mkdir()
    egs_home.mkdir()

    backend = EgsnrcReferenceBackend(
        beamnrc_executable=beamnrc,
        dosxyznrc_executable=dosxyz,
        hen_house=hen_house,
        egs_home=egs_home,
    )

    assert backend.validate_environment() == ()
    assert backend.resolved_beamnrc() == beamnrc
    assert backend.resolved_dosxyznrc() == dosxyz


def test_environment_validation_reports_missing_components(tmp_path: Path) -> None:
    backend = EgsnrcReferenceBackend(
        beamnrc_executable=tmp_path / "missing-beamnrc",
        dosxyznrc_executable=tmp_path / "missing-dosxyz",
        hen_house=tmp_path / "missing-hen-house",
        egs_home=tmp_path / "missing-egs-home",
    )

    errors = backend.validate_environment()

    assert any("BEAMnrc executable not found" in item for item in errors)
    assert any("DOSXYZnrc executable not found" in item for item in errors)
    assert any("HEN_HOUSE directory does not exist" in item for item in errors)
    assert any("EGS_HOME directory does not exist" in item for item in errors)
