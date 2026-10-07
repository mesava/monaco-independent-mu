from __future__ import annotations

import numpy as np
import pytest

from indep_mu.dicom.rtplan import (
    Beam,
    BeamLimitingDeviceDefinition,
    BeamLimitingDeviceState,
    ControlPoint,
)
from indep_mu.montecarlo.source21 import (
    Source21Angles,
    build_source21_control_points,
    patient_point_to_egsphant_cm,
    render_source21_control_points,
    validate_synchronized_mu_indices,
)
from indep_mu.montecarlo.syncmlce import (
    SourceFocusedSyncMlceMapper,
    build_syncmlce_sequence,
)
from indep_mu.patient_model.egsphant import EgsphantGeometry


class DummyOrientationMapper:
    def map_angles(
        self,
        *,
        gantry_deg: float,
        collimator_deg: float,
        patient_support_deg: float,
    ) -> Source21Angles:
        return Source21Angles(
            theta_deg=gantry_deg,
            phi_deg=patient_support_deg,
            phicol_deg=collimator_deg,
        )


def _geometry() -> EgsphantGeometry:
    return EgsphantGeometry(
        x_bound_cm=np.asarray([-1.0, 0.0, 1.0]),
        y_bound_cm=np.asarray([-2.0, 0.0, 2.0]),
        z_bound_cm=np.asarray([-3.0, 0.0, 3.0]),
        x_direction_patient=(1.0, 0.0, 0.0),
        y_direction_patient=(0.0, 1.0, 0.0),
        z_direction_patient=(0.0, 0.0, 1.0),
    )


def _cp(index: int, cmw: float, mlc: tuple[float, ...]) -> ControlPoint:
    return ControlPoint(
        index=index,
        cumulative_meterset_weight=cmw,
        gantry_angle_deg=10.0 + index,
        gantry_rotation_direction="CC",
        collimator_angle_deg=20.0,
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
        isocenter_position_mm=(10.0, 20.0, 30.0),
        device_positions=(
            BeamLimitingDeviceState("MLCX", mlc),
        ),
    )


def _beam(cmw=(0.0, 0.5, 1.0)) -> Beam:
    cps = tuple(
        _cp(
            index,
            value,
            (-10.0 + index, -8.0 + index, 8.0 - index, 10.0 - index),
        )
        for index, value in enumerate(cmw)
    )
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
        ),
        control_points=cps,
    )


def test_patient_point_uses_egsphant_axis_basis() -> None:
    geometry = EgsphantGeometry(
        x_bound_cm=np.asarray([0.0, 1.0]),
        y_bound_cm=np.asarray([0.0, 1.0]),
        z_bound_cm=np.asarray([0.0, 1.0]),
        x_direction_patient=(0.0, 1.0, 0.0),
        y_direction_patient=(1.0, 0.0, 0.0),
        z_direction_patient=(0.0, 0.0, -1.0),
    )

    point = patient_point_to_egsphant_cm((10.0, 20.0, 30.0), geometry)
    assert point == pytest.approx((2.0, 1.0, -3.0))


def test_source21_control_points_and_rendering() -> None:
    source = build_source21_control_points(
        _beam(),
        phantom_geometry=_geometry(),
        orientation_mapper=DummyOrientationMapper(),
        dsource_cm=50.0,
    )

    assert len(source.points) == 3
    assert source.points[0].mu_index == pytest.approx(0.0)
    assert source.points[-1].mu_index == pytest.approx(1.0)
    assert source.points[0].xiso_cm == pytest.approx(1.0)
    assert source.points[0].yiso_cm == pytest.approx(2.0)
    assert source.points[0].ziso_cm == pytest.approx(3.0)

    rendered = render_source21_control_points(source)
    assert "    50.00000000" in rendered


def test_source21_and_syncmlce_share_mu_grid() -> None:
    beam = _beam()
    source = build_source21_control_points(
        beam,
        phantom_geometry=_geometry(),
        orientation_mapper=DummyOrientationMapper(),
        dsource_cm=50.0,
    )
    mlc = build_syncmlce_sequence(
        beam,
        mapper=SourceFocusedSyncMlceMapper(
            sad_mm=1000.0,
            zmin_cm=40.0,
            negative_bank=1,
        ),
    )

    validate_synchronized_mu_indices(source, mlc=mlc)


def test_mu_grid_mismatch_is_rejected() -> None:
    source = build_source21_control_points(
        _beam(),
        phantom_geometry=_geometry(),
        orientation_mapper=DummyOrientationMapper(),
        dsource_cm=50.0,
    )
    mlc = build_syncmlce_sequence(
        _beam(cmw=(0.0, 0.4, 1.0)),
        mapper=SourceFocusedSyncMlceMapper(
            sad_mm=1000.0,
            zmin_cm=40.0,
            negative_bank=1,
        ),
    )

    with pytest.raises(ValueError, match="not synchronized"):
        validate_synchronized_mu_indices(source, mlc=mlc)
