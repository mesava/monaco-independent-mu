import numpy as np
import pytest

from indep_mu.beam_model.agility_rounded_tip import (
    RoundedLeafTipTangentGeometry,
    agility_literature_tangent_geometry,
    leaf_bank_rotation_shift_mm,
)


@pytest.mark.parametrize("side", ["negative", "positive"])
@pytest.mark.parametrize(
    "x_iso_cm",
    [-20.0, -5.0, 0.0, 3.25, 20.0],
)
def test_tangent_mapping_round_trip(side: str, x_iso_cm: float) -> None:
    geometry = agility_literature_tangent_geometry()

    cylinder = geometry.projected_edge_to_cylinder_origin_cm(
        x_iso_cm,
        opening_side=side,
    )
    restored = geometry.cylinder_origin_to_projected_edge_cm(
        cylinder,
        opening_side=side,
    )

    assert float(restored) == pytest.approx(x_iso_cm, abs=1e-12)


@pytest.mark.parametrize("side", ["negative", "positive"])
@pytest.mark.parametrize("x_iso_cm", [-15.0, -2.0, 0.0, 8.0, 15.0])
def test_source_ray_is_exactly_tangent(side: str, x_iso_cm: float) -> None:
    geometry = agility_literature_tangent_geometry()
    cylinder = geometry.projected_edge_to_cylinder_origin_cm(
        x_iso_cm,
        opening_side=side,
    )

    distance = geometry.tangent_distance_cm(x_iso_cm, cylinder)

    assert float(distance) == pytest.approx(geometry.radius_cm, abs=1e-12)


def test_zero_projected_edge_has_expected_cylinder_centres() -> None:
    geometry = agility_literature_tangent_geometry()

    negative = geometry.projected_edge_to_cylinder_origin_cm(
        0.0,
        opening_side="negative",
    )
    positive = geometry.projected_edge_to_cylinder_origin_cm(
        0.0,
        opening_side="positive",
    )

    assert float(negative) == pytest.approx(-17.0)
    assert float(positive) == pytest.approx(17.0)


def test_vectorized_round_trip() -> None:
    geometry = RoundedLeafTipTangentGeometry(
        sad_cm=100.0,
        radius_cm=17.0,
        cylinder_axis_z_cm=34.93,
    )
    projected = np.linspace(-20.0, 20.0, 17)

    for side in ("negative", "positive"):
        cylinder = geometry.projected_edge_to_cylinder_origin_cm(
            projected,
            opening_side=side,
        )
        restored = geometry.cylinder_origin_to_projected_edge_cm(
            cylinder,
            opening_side=side,
        )
        np.testing.assert_allclose(restored, projected, rtol=0.0, atol=1e-11)


def test_invalid_geometry_is_rejected() -> None:
    with pytest.raises(ValueError, match="farther from the source"):
        RoundedLeafTipTangentGeometry(
            sad_cm=100.0,
            radius_cm=17.0,
            cylinder_axis_z_cm=10.0,
        )



def test_published_agility_lbrot_shift() -> None:
    shift_mm = leaf_bank_rotation_shift_mm(
        leaf_thickness_mm=90.0,
        lbrot_rad=0.009,
    )

    # Equation gives about 0.405 mm; the publication reports 0.41 mm from
    # the formula and about 0.42 mm from measurement/MC.
    assert shift_mm == pytest.approx(0.405, abs=2e-5)
    assert shift_mm == pytest.approx(0.41, abs=0.01)



@pytest.mark.parametrize("side", ["negative", "positive"])
def test_agility_clinical_edge_range_tangent_stays_within_leaf_slab(
    side: str,
) -> None:
    geometry = agility_literature_tangent_geometry()
    projected = np.linspace(-20.0, 20.0, 81)

    cylinder = geometry.projected_edge_to_cylinder_origin_cm(
        projected,
        opening_side=side,
    )
    inside = geometry.tangent_within_leaf_slab(
        projected,
        cylinder,
        zmin_cm=31.18,
        zmax_cm=40.18,
    )

    assert np.all(inside)


def test_extreme_projected_edge_can_leave_rounded_tip_slab() -> None:
    geometry = agility_literature_tangent_geometry()

    projected = np.asarray([40.0])
    cylinder = geometry.projected_edge_to_cylinder_origin_cm(
        projected,
        opening_side="negative",
    )

    inside = geometry.tangent_within_leaf_slab(
        projected,
        cylinder,
        zmin_cm=31.18,
        zmax_cm=40.18,
    )

    assert not bool(inside[0])
