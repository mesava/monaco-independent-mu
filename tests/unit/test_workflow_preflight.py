import numpy as np

from indep_mu.dicom.ct import CtGeometry, CtSeries
from indep_mu.dicom.rtplan import (
    Beam,
    ControlPoint,
    RtPlan,
)
from indep_mu.workflow.preflight import run_transport_preflight


def _ct(*, patient_position: str | None = "HFS", kvp: float | None = 120.0) -> CtSeries:
    geometry = CtGeometry(
        rows=2,
        columns=2,
        pixel_spacing_mm=(1.0, 1.0),
        slice_spacing_mm=2.5,
        image_orientation_patient=(1.0, 0.0, 0.0, 0.0, 1.0, 0.0),
        first_image_position_patient_mm=(0.0, 0.0, 0.0),
        last_image_position_patient_mm=(0.0, 0.0, 2.5),
        patient_position=patient_position,
        frame_of_reference_uid="1.2.3",
        series_instance_uid="1.2.3.4",
        kvp=kvp,
    )
    return CtSeries(
        hu=np.zeros((2, 2, 2), dtype=np.float32),
        geometry=geometry,
        slice_positions_mm=np.asarray([0.0, 2.5]),
        image_positions_patient_mm=np.asarray(
            [[0.0, 0.0, 0.0], [0.0, 0.0, 2.5]]
        ),
        sop_instance_uids=("1", "2"),
    )


def _cp(index: int, gantry: float, couch: float) -> ControlPoint:
    return ControlPoint(
        index=index,
        cumulative_meterset_weight=float(index),
        gantry_angle_deg=gantry,
        gantry_rotation_direction="NONE",
        collimator_angle_deg=0.0,
        collimator_rotation_direction="NONE",
        patient_support_angle_deg=couch,
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
        device_positions=(),
    )


def _plan(*, gantry: float = 0.0, couch: float = 0.0) -> RtPlan:
    beam = Beam(
        number=1,
        name="TEST",
        description=None,
        beam_type="STATIC",
        radiation_type="PHOTON",
        treatment_delivery_type="TREATMENT",
        treatment_machine_name="VERSA",
        source_axis_distance_mm=1000.0,
        final_cumulative_meterset_weight=1.0,
        beam_meterset_mu=100.0,
        fluence_mode="STANDARD",
        fluence_mode_id=None,
        number_of_wedges=0,
        number_of_compensators=0,
        number_of_boli=0,
        number_of_blocks=0,
        device_definitions=(),
        control_points=(
            _cp(0, gantry, couch),
            ControlPoint(
                **{
                    **_cp(1, gantry, couch).__dict__,
                    "cumulative_meterset_weight": 1.0,
                }
            ),
        ),
    )
    return RtPlan(
        sop_instance_uid="9.8.7",
        plan_label="TEST",
        plan_name="TEST",
        fraction_group_number=1,
        number_of_fractions_planned=1,
        beams=(beam,),
    )


def test_transport_preflight_passes_hfs_120kv() -> None:
    report = run_transport_preflight(_ct(), _plan())
    assert report.passed


def test_transport_preflight_rejects_non_hfs() -> None:
    report = run_transport_preflight(_ct(patient_position="FFS"), _plan())
    assert not report.passed
    assert any(item.code == "PATIENT_POSITION_UNSUPPORTED" for item in report.errors)


def test_transport_preflight_rejects_kvp_mismatch() -> None:
    report = run_transport_preflight(_ct(kvp=100.0), _plan())
    assert not report.passed
    assert any(item.code == "CT_CALIBRATION_KVP_MISMATCH" for item in report.errors)


def test_transport_preflight_rejects_source21_singularity() -> None:
    report = run_transport_preflight(
        _ct(),
        _plan(gantry=90.0, couch=90.0),
    )
    assert not report.passed
    assert any(
        item.code == "SOURCE21_ORIENTATION_UNSUPPORTED"
        for item in report.errors
    )
