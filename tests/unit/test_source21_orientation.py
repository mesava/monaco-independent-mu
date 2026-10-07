import pytest

from indep_mu.montecarlo.source21_orientation import (
    ZhanHfsSource21OrientationMapper,
)
from indep_mu.montecarlo.source21_reference import (
    coplanar_hfs_reference_angles,
)


@pytest.mark.parametrize(
    ("gantry", "collimator"),
    [
        (0.0, 0.0),
        (90.0, 0.0),
        (180.0, 0.0),
        (270.0, 0.0),
        (30.0, 15.0),
        (240.0, 90.0),
    ],
)
def test_zhan_mapper_matches_coplanar_reference(
    gantry: float,
    collimator: float,
) -> None:
    mapper = ZhanHfsSource21OrientationMapper()
    actual = mapper.map_angles(
        gantry_deg=gantry,
        collimator_deg=collimator,
        patient_support_deg=0.0,
    )
    expected = coplanar_hfs_reference_angles(
        gantry_deg=gantry,
        collimator_deg=collimator,
        patient_support_deg=0.0,
    )

    assert actual.theta_deg == pytest.approx(expected.theta_deg, abs=1e-12)
    assert actual.phi_deg == pytest.approx(expected.phi_deg, abs=1e-12)
    assert actual.phicol_deg == pytest.approx(expected.phicol_deg, abs=1e-12)


@pytest.mark.parametrize(
    ("gantry", "couch", "collimator", "theta", "phi", "phicol"),
    [
        (
            30.0,
            20.0,
            15.0,
            99.84655193983407,
            298.4812382813395,
            272.49524075699975,
        ),
        (
            120.0,
            15.0,
            25.0,
            102.95253964222236,
            30.86747779067436,
            237.36925978756994,
        ),
        (
            240.0,
            330.0,
            90.0,
            115.6589062732553,
            146.30993247402017,
            196.10211375198605,
        ),
    ],
)
def test_zhan_non_coplanar_regression_vectors(
    gantry: float,
    couch: float,
    collimator: float,
    theta: float,
    phi: float,
    phicol: float,
) -> None:
    mapper = ZhanHfsSource21OrientationMapper()
    actual = mapper.map_angles(
        gantry_deg=gantry,
        collimator_deg=collimator,
        patient_support_deg=couch,
    )

    assert actual.theta_deg == pytest.approx(theta, abs=1e-10)
    assert actual.phi_deg == pytest.approx(phi, abs=1e-10)
    assert actual.phicol_deg == pytest.approx(phicol, abs=1e-10)


@pytest.mark.parametrize(
    ("gantry", "couch"),
    [
        (90.0, 90.0),
        (90.0, 270.0),
        (270.0, 90.0),
        (270.0, 270.0),
    ],
)
def test_exact_singularities_are_rejected(
    gantry: float,
    couch: float,
) -> None:
    mapper = ZhanHfsSource21OrientationMapper()

    with pytest.raises(ValueError, match="singular"):
        mapper.map_angles(
            gantry_deg=gantry,
            collimator_deg=0.0,
            patient_support_deg=couch,
        )


def test_non_hfs_is_not_silently_accepted() -> None:
    with pytest.raises(ValueError, match="only for HFS"):
        ZhanHfsSource21OrientationMapper(patient_position="FFS")
