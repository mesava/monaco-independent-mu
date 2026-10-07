import pytest

from indep_mu.montecarlo.source21_reference import (
    coplanar_hfs_reference_angles,
)


@pytest.mark.parametrize(
    ("gantry", "expected_phi"),
    [
        (0.0, 270.0),
        (90.0, 0.0),
        (180.0, 90.0),
        (270.0, 180.0),
    ],
)
def test_coplanar_cardinal_gantry_reference(
    gantry: float,
    expected_phi: float,
) -> None:
    result = coplanar_hfs_reference_angles(
        gantry_deg=gantry,
        collimator_deg=0.0,
        patient_support_deg=0.0,
    )

    assert result.theta_deg == pytest.approx(90.0)
    assert result.phi_deg == pytest.approx(expected_phi)
    assert result.phicol_deg == pytest.approx(270.0)


def test_coplanar_collimator_reference() -> None:
    result = coplanar_hfs_reference_angles(
        gantry_deg=0.0,
        collimator_deg=90.0,
        patient_support_deg=0.0,
    )
    assert result.phicol_deg == pytest.approx(180.0)


def test_nonzero_couch_is_rejected() -> None:
    with pytest.raises(ValueError, match="couch angle 0"):
        coplanar_hfs_reference_angles(
            gantry_deg=0.0,
            collimator_deg=0.0,
            patient_support_deg=15.0,
        )
