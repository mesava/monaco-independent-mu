from pathlib import Path

import pytest

from indep_mu.beam_model.tps_reference import load_tps_reference_yaml


def test_loads_monaco_reference_values(tmp_path: Path) -> None:
    path = tmp_path / "reference.yaml"
    path.write_text(
        """
machine_model_id: VERSA_HD
tps_name: Elekta Monaco
tps_version: "6.1.4"
reference_geometry:
  field: {x_mm: 100, y_mm: 100}
  ssd_mm: 900
  depth_mm: 100
  delivered_mu: 100
energies:
  - {energy_model_id: 6MV, dose_gy: 0.995}
  - {energy_model_id: 6FFF, dose_gy: 0.994}
  - {energy_model_id: 10MV, dose_gy: 1.000}
  - {energy_model_id: 10FFF, dose_gy: 1.001}
""",
        encoding="utf-8",
    )

    data = load_tps_reference_yaml(path)

    assert data.geometry.delivered_mu == pytest.approx(100.0)
    assert data.dose_gy_per_100mu("6MV") == pytest.approx(0.995)
    assert data.dose_gy_per_100mu("6FFF") == pytest.approx(0.994)
    assert data.dose_gy_per_100mu("10MV") == pytest.approx(1.000)
    assert data.dose_gy_per_100mu("10FFF") == pytest.approx(1.001)
    assert data.dose_gy_per_mu("6MV") == pytest.approx(0.00995)


def test_requires_explicit_delivered_mu(tmp_path: Path) -> None:
    path = tmp_path / "bad.yaml"
    path.write_text(
        """
machine_model_id: VERSA_HD
tps_name: Elekta Monaco
reference_geometry:
  field: {x_mm: 100, y_mm: 100}
  ssd_mm: 900
  depth_mm: 100
energies:
  - {energy_model_id: 6MV, dose_gy: 0.995}
""",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="delivered_mu"):
        load_tps_reference_yaml(path)
