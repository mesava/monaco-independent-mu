from pathlib import Path

import pytest

from indep_mu.beam_model.commissioning import (
    FieldSize,
    ReferenceCalibration,
    load_commissioning_yaml,
)


def test_reference_calibration_never_assumes_mu() -> None:
    calibration = ReferenceCalibration(
        field=FieldSize(100.0, 100.0),
        ssd_mm=900.0,
        depth_mm=100.0,
        delivered_mu=100.0,
        measured_dose_gy=0.995,
    )

    assert calibration.dose_gy_per_mu == pytest.approx(0.00995)
    assert calibration.dose_gy_per_100mu == pytest.approx(0.995)


def test_yaml_requires_explicit_delivered_mu(tmp_path: Path) -> None:
    path = tmp_path / "commissioning.yaml"
    path.write_text(
        """
machine_model_id: VERSA_HD
energies:
  - energy_model_id: 6MV
    calibration:
      field: {x_mm: 100, y_mm: 100}
      ssd_mm: 900
      depth_mm: 100
      measured_dose_gy: 0.995
""",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="delivered_mu"):
        load_commissioning_yaml(path)


def test_yaml_loads_output_factor_and_curves(tmp_path: Path) -> None:
    path = tmp_path / "commissioning.yaml"
    path.write_text(
        """
machine_model_id: VERSA_HD
energies:
  - energy_model_id: 6MV
    calibration:
      field: {x_mm: 100, y_mm: 100}
      ssd_mm: 900
      depth_mm: 100
      delivered_mu: 100
      measured_dose_gy: 0.995
    output_factors:
      - field: {x_mm: 100, y_mm: 100}
        output_factor: 1.0
    pdd_curves:
      - field: {x_mm: 100, y_mm: 100}
        ssd_mm: 1000
        depth_mm: [0, 10, 20]
        relative_dose: [0.5, 1.0, 0.9]
    profiles:
      - field: {x_mm: 100, y_mm: 100}
        ssd_mm: 1000
        depth_mm: 100
        axis: CROSSPLANE
        position_mm: [-50, 0, 50]
        relative_dose: [0.5, 1.0, 0.5]
""",
        encoding="utf-8",
    )

    dataset = load_commissioning_yaml(path)
    energy = dataset.energy("6MV")

    assert energy.calibration.dose_gy_per_100mu == pytest.approx(0.995)
    assert len(energy.output_factors) == 1
    assert len(energy.pdd_curves) == 1
    assert len(energy.profiles) == 1
    energy.validate_reference_output_factor()
