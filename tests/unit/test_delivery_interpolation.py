from pathlib import Path

import pytest

from indep_mu.dicom.rtplan import load_rtplan
from indep_mu.montecarlo.delivery import (
    DeliveryInterpolationPolicy,
    directed_rotation_delta_deg,
    sample_beam_segment,
)
from test_rtplan import _write_plan


POLICY = DeliveryInterpolationPolicy(
    policy_id="synthetic-linear",
    evidence="Unit-test-only linear interpolation policy.",
)


def test_dicom_rotation_examples() -> None:
    assert directed_rotation_delta_deg(5, 5, "NONE") == pytest.approx(0.0)
    assert directed_rotation_delta_deg(5, 5, "CW") == pytest.approx(-360.0)
    assert directed_rotation_delta_deg(170, 160, "CC") == pytest.approx(350.0)


def test_none_direction_rejects_changed_angle() -> None:
    with pytest.raises(ValueError, match="NONE"):
        directed_rotation_delta_deg(10, 20, "NONE")


def test_sample_segment_interpolates_cmw_mu_and_mlc(tmp_path: Path) -> None:
    path = tmp_path / "plan.dcm"
    _write_plan(path)
    beam = load_rtplan(path).beam(1)

    state = sample_beam_segment(
        beam,
        0,
        0.5,
        policy=POLICY,
    )

    assert state.cumulative_meterset_weight == pytest.approx(0.2)
    assert state.mu_from_beam_start == pytest.approx(100.0)
    assert state.positions("ASYMY") == pytest.approx((-50.0, 50.0))
    assert state.positions("MLCX") == pytest.approx(
        (-9.5, -7.5, 7.5, 9.5)
    )
    # CP0 encodes 180 -> 150 with CC. DICOM machine-rotation semantics mean
    # increasing IEC angle, i.e. a 330-degree CC path; the midpoint is 345 deg.
    assert state.gantry_angle_deg == pytest.approx(345.0)
