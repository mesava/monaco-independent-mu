from __future__ import annotations

import numpy as np
import pytest

from indep_mu.beam_model.agility_rounded_tip import RoundedLeafTipTangentGeometry
from indep_mu.dicom.rtplan import (
    Beam,
    BeamLimitingDeviceDefinition,
    BeamLimitingDeviceState,
    ControlPoint,
)
from indep_mu.montecarlo.syncmlce import (
    ResearchRoundedTipSyncMlceMapper,
    SourceFocusedSyncMlceMapper,
    build_syncmlce_sequence,
    render_syncmlce_sequence,
)
from indep_mu.montecarlo.syncjaws import (
    FocusedJawPairGeometry,
    build_syncjaws_sequence,
    render_syncjaws_sequence,
)


def _cp(
    index: int,
    cmw: float,
    *,
    mlc: tuple[float, ...],
    jaw: tuple[float, ...],
) -> ControlPoint:
    return ControlPoint(
        index=index,
        cumulative_meterset_weight=cmw,
        gantry_angle_deg=180.0,
        gantry_rotation_direction="NONE",
        collimator_angle_deg=0.0,
        collimator_rotation_direction="NONE",
        patient_support_angle_deg=0.0,
        patient_support_rotation_direction="NONE",
        table_top_eccentric_angle_deg=None,
        table_top_eccentric_rotation_direction=None,
        table_top_vertical_position_mm=None,
        table_top_longitudinal_position_mm=None,
        table_top_lateral_position_mm=None,
        nominal_beam_energy_mv=6.0,
        dose_rate_set_mu_min=600.0,
        source_to_surface_distance_mm=900.0,
        isocenter_position_mm=(0.0, 0.0, 0.0),
        device_positions=(
            BeamLimitingDeviceState("MLCX", mlc),
            BeamLimitingDeviceState("ASYMY", jaw),
        ),
    )


def _beam(control_points: tuple[ControlPoint, ...]) -> Beam:
    return Beam(
        number=1,
        name="ARC1",
        description=None,
        beam_type="DYNAMIC",
        radiation_type="PHOTON",
        treatment_delivery_type="TREATMENT",
        treatment_machine_name="VERSA",
        source_axis_distance_mm=1000.0,
        final_cumulative_meterset_weight=1.0,
        beam_meterset_mu=200.0,
        fluence_mode="STANDARD",
        fluence_mode_id=None,
        number_of_wedges=0,
        number_of_compensators=0,
        number_of_boli=0,
        number_of_blocks=0,
        device_definitions=(
            BeamLimitingDeviceDefinition(
                device_type="MLCX",
                number_of_leaf_jaw_pairs=2,
                leaf_position_boundaries_mm=(-5.0, 0.0, 5.0),
                source_to_device_distance_mm=400.0,
            ),
            BeamLimitingDeviceDefinition(
                device_type="ASYMY",
                number_of_leaf_jaw_pairs=1,
                leaf_position_boundaries_mm=None,
                source_to_device_distance_mm=450.0,
            ),
        ),
        control_points=control_points,
    )


def test_syncmlce_focused_projection_and_rendering() -> None:
    beam = _beam(
        (
            _cp(
                0,
                0.0,
                mlc=(-10.0, -8.0, 8.0, 10.0),
                jaw=(-50.0, 50.0),
            ),
            _cp(
                1,
                1.0,
                mlc=(-5.0, -4.0, 5.0, 6.0),
                jaw=(-40.0, 40.0),
            ),
        )
    )

    sequence = build_syncmlce_sequence(
        beam,
        mapper=SourceFocusedSyncMlceMapper(
            sad_mm=1000.0,
            zmin_cm=40.0,
            negative_bank=1,
        ),
    )

    np.testing.assert_allclose(
        sequence.points[0].opening.negative_cm,
        [-0.4, -0.32],
    )
    np.testing.assert_allclose(
        sequence.points[0].opening.positive_cm,
        [0.32, 0.4],
    )

    rendered = render_syncmlce_sequence(sequence)
    assert "ARC1" in rendered
    assert f"{2:10d}" in rendered
    assert "     0.00000000" in rendered
    assert "     1.00000000" in rendered


def test_equal_cmw_keeps_last_zero_weight_state() -> None:
    beam = _beam(
        (
            _cp(
                0,
                0.0,
                mlc=(-20.0, -20.0, 20.0, 20.0),
                jaw=(-50.0, 50.0),
            ),
            _cp(
                1,
                0.0,
                mlc=(-10.0, -8.0, 8.0, 10.0),
                jaw=(-50.0, 50.0),
            ),
            _cp(
                2,
                1.0,
                mlc=(-5.0, -4.0, 5.0, 6.0),
                jaw=(-40.0, 40.0),
            ),
        )
    )

    sequence = build_syncmlce_sequence(
        beam,
        mapper=SourceFocusedSyncMlceMapper(
            sad_mm=1000.0,
            zmin_cm=40.0,
            negative_bank=1,
        ),
    )

    assert len(sequence.points) == 2
    np.testing.assert_allclose(
        sequence.points[0].opening.negative_cm,
        [-0.4, -0.32],
    )


def test_syncjaws_projects_front_and_back_surfaces() -> None:
    beam = _beam(
        (
            _cp(
                0,
                0.0,
                mlc=(-10.0, -8.0, 8.0, 10.0),
                jaw=(-50.0, 50.0),
            ),
            _cp(
                1,
                1.0,
                mlc=(-5.0, -4.0, 5.0, 6.0),
                jaw=(-40.0, 40.0),
            ),
        )
    )

    sequence = build_syncjaws_sequence(
        beam,
        geometries=(
            FocusedJawPairGeometry(
                device_type="ASYMY",
                zmin_cm=40.0,
                zmax_cm=50.0,
                negative_bank=1,
            ),
        ),
    )

    first = sequence.points[0].pairs[0]
    assert first.negative_front_cm == pytest.approx(-2.0)
    assert first.positive_front_cm == pytest.approx(2.0)
    assert first.negative_back_cm == pytest.approx(-2.5)
    assert first.positive_back_cm == pytest.approx(2.5)

    rendered = render_syncjaws_sequence(sequence)
    assert f"{2:10d}" in rendered
    assert "     0.00000000" in rendered
    assert "     1.00000000" in rendered



def test_research_rounded_tip_mapper_recovers_projected_edges() -> None:
    beam = _beam(
        (
            _cp(
                0,
                0.0,
                mlc=(-50.0, -30.0, 40.0, 60.0),
                jaw=(-50.0, 50.0),
            ),
            _cp(
                1,
                1.0,
                mlc=(-40.0, -20.0, 30.0, 50.0),
                jaw=(-40.0, 40.0),
            ),
        )
    )

    mapper = ResearchRoundedTipSyncMlceMapper(
        sad_mm=1000.0,
        radius_cm=17.0,
        cylinder_axis_z_cm=34.93,
        negative_bank=1,
    )
    sequence = build_syncmlce_sequence(beam, mapper=mapper)

    geometry = RoundedLeafTipTangentGeometry(
        sad_cm=100.0,
        radius_cm=17.0,
        cylinder_axis_z_cm=34.93,
    )

    first = sequence.points[0].opening
    negative_edge = geometry.cylinder_origin_to_projected_edge_cm(
        first.negative_cm,
        opening_side="negative",
    )
    positive_edge = geometry.cylinder_origin_to_projected_edge_cm(
        first.positive_cm,
        opening_side="positive",
    )

    np.testing.assert_allclose(negative_edge, [-5.0, -3.0], atol=1e-11)
    np.testing.assert_allclose(positive_edge, [4.0, 6.0], atol=1e-11)


def test_research_rounded_mapper_bank_identity_is_explicit() -> None:
    from indep_mu.beam_model.iec_coordinates import DicomBankPositions

    banks = DicomBankPositions(
        bank1_mm=np.asarray([30.0]),
        bank2_mm=np.asarray([-40.0]),
    )
    mapper = ResearchRoundedTipSyncMlceMapper(
        sad_mm=1000.0,
        radius_cm=17.0,
        cylinder_axis_z_cm=34.93,
        negative_bank=2,
    )

    opening = mapper.map_banks(banks)

    geometry = RoundedLeafTipTangentGeometry(
        sad_cm=100.0,
        radius_cm=17.0,
        cylinder_axis_z_cm=34.93,
    )
    restored_negative = geometry.cylinder_origin_to_projected_edge_cm(
        opening.negative_cm,
        opening_side="negative",
    )
    restored_positive = geometry.cylinder_origin_to_projected_edge_cm(
        opening.positive_cm,
        opening_side="positive",
    )

    np.testing.assert_allclose(restored_negative, [-4.0], atol=1e-11)
    np.testing.assert_allclose(restored_positive, [3.0], atol=1e-11)



def test_research_rounded_mapper_applies_explicit_projected_shift() -> None:
    from indep_mu.beam_model.iec_coordinates import DicomBankPositions

    banks = DicomBankPositions(
        bank1_mm=np.asarray([-50.0]),
        bank2_mm=np.asarray([50.0]),
    )
    mapper = ResearchRoundedTipSyncMlceMapper(
        sad_mm=1000.0,
        radius_cm=17.0,
        cylinder_axis_z_cm=34.93,
        negative_bank=1,
        projected_edge_shift_mm=0.405,
    )

    opening = mapper.map_banks(banks)
    geometry = RoundedLeafTipTangentGeometry(
        sad_cm=100.0,
        radius_cm=17.0,
        cylinder_axis_z_cm=34.93,
    )

    negative_edge_cm = geometry.cylinder_origin_to_projected_edge_cm(
        opening.negative_cm,
        opening_side="negative",
    )
    positive_edge_cm = geometry.cylinder_origin_to_projected_edge_cm(
        opening.positive_cm,
        opening_side="positive",
    )

    np.testing.assert_allclose(
        negative_edge_cm,
        [(-50.0 + 0.405) / 10.0],
        atol=1e-11,
    )
    np.testing.assert_allclose(
        positive_edge_cm,
        [(50.0 + 0.405) / 10.0],
        atol=1e-11,
    )



def test_research_rounded_mapper_checks_physical_leaf_slab() -> None:
    from indep_mu.beam_model.iec_coordinates import DicomBankPositions

    mapper = ResearchRoundedTipSyncMlceMapper(
        sad_mm=1000.0,
        radius_cm=17.0,
        cylinder_axis_z_cm=34.93,
        negative_bank=1,
        zmin_cm=31.18,
        zmax_cm=40.18,
    )

    normal = DicomBankPositions(
        bank1_mm=np.asarray([-200.0]),
        bank2_mm=np.asarray([200.0]),
    )
    opening = mapper.map_banks(normal)
    assert opening.negative_cm.shape == (1,)
    assert opening.positive_cm.shape == (1,)

    extreme = DicomBankPositions(
        bank1_mm=np.asarray([400.0]),
        bank2_mm=np.asarray([200.0]),
    )
    with pytest.raises(ValueError, match="tangent outside"):
        mapper.map_banks(extreme)
