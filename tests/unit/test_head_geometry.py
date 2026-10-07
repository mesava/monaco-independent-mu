from pathlib import Path

import yaml

from indep_mu.beam_model.head_geometry import (
    load_head_geometry_yaml,
    validate_head_geometry,
)


def test_incomplete_head_geometry_is_blocked(tmp_path: Path) -> None:
    path = tmp_path / "head.yaml"
    path.write_text(
        """
machine_model_id: VERSA_HD
mlc:
  end_type: rounded
  zmin_cm: null
  zmax_cm: null
  leaf_radius_cm: null
  cylinder_axis_z_cm: null
  leaf_bank_tilt_rad: null
  negative_bank: null
jaws:
  - device_type: ASYMY
    zmin_cm: null
    zmax_cm: null
    negative_bank: null
evidence: []
""",
        encoding="utf-8",
    )

    geometry = load_head_geometry_yaml(path)
    report = validate_head_geometry(geometry)

    assert report.errors
    assert not report.ready_for_source_focused_jaws
    assert not report.ready_for_rounded_agility_syncmlce
    assert any(
        item.code == "MLC_ROUNDED_MAPPING_UNVALIDATED"
        for item in report.errors
    )


def test_complete_numbers_still_do_not_enable_unvalidated_rounded_mapping(
    tmp_path: Path,
) -> None:
    path = tmp_path / "head.yaml"
    raw = {
        "machine_model_id": "VERSA_HD",
        "mlc": {
            "end_type": "rounded",
            "zmin_cm": 35.0,
            "zmax_cm": 42.0,
            "leaf_radius_cm": 10.0,
            "cylinder_axis_z_cm": 38.0,
            "leaf_bank_tilt_rad": 0.0,
            "negative_bank": 1,
        },
        "jaws": [
            {
                "device_type": "ASYMY",
                "zmin_cm": 30.0,
                "zmax_cm": 34.0,
                "negative_bank": 1,
            }
        ],
        "evidence": ["synthetic test values only"],
    }
    path.write_text(yaml.safe_dump(raw), encoding="utf-8")

    report = validate_head_geometry(load_head_geometry_yaml(path))

    assert report.ready_for_source_focused_jaws
    assert not report.ready_for_rounded_agility_syncmlce
    assert [item.code for item in report.errors] == [
        "MLC_ROUNDED_MAPPING_UNVALIDATED"
    ]
