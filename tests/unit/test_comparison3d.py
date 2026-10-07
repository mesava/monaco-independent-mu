import numpy as np
import pytest

from indep_mu.comparison3d import (
    PatientMcDoseGridGy,
    compare_on_rtdose_grid,
    dvh_metrics,
)
from indep_mu.comparison3d.dose_quantity import DoseQuantity
from indep_mu.comparison3d.dvh import StructureDoseSamples
from indep_mu.dicom.rtdose import DoseGeometry, RtDose
from indep_mu.montecarlo.threedose import ThreeDDose
from indep_mu.patient_model.egsphant import EgsphantGeometry


def _mc_grid() -> PatientMcDoseGridGy:
    boundaries = np.asarray([-1.5, -0.5, 0.5, 1.5], dtype=np.float64)
    centres = 0.5 * (boundaries[:-1] + boundaries[1:])
    zz, yy, xx = np.meshgrid(centres, centres, centres, indexing="ij")
    dose = 2.0 + 0.1 * xx + 0.2 * yy + 0.3 * zz

    return PatientMcDoseGridGy(
        dose=ThreeDDose(
            x_bound_cm=boundaries,
            y_bound_cm=boundaries,
            z_bound_cm=boundaries,
            dose=dose,
            relative_uncertainty=np.full(dose.shape, 0.01),
        ),
        phantom_geometry=EgsphantGeometry(
            x_bound_cm=boundaries,
            y_bound_cm=boundaries,
            z_bound_cm=boundaries,
            x_direction_patient=(1.0, 0.0, 0.0),
            y_direction_patient=(0.0, 1.0, 0.0),
            z_direction_patient=(0.0, 0.0, 1.0),
        ),
        frame_of_reference_uid="1.2.3",
        dose_quantity=DoseQuantity.DOSE_TO_MEDIUM,
    )


def test_mc_patient_coordinate_sampling() -> None:
    grid = _mc_grid()

    points_mm = np.asarray([[2.5, -5.0, 7.5]])
    value = grid.sample_patient_points_gy(points_mm)[0]

    expected = 2.0 + 0.1 * 0.25 + 0.2 * (-0.5) + 0.3 * 0.75
    assert value == pytest.approx(expected)


def test_dvh_percentile_semantics() -> None:
    samples = StructureDoseSamples(
        dose_gy=np.arange(1.0, 101.0),
        voxel_volume_cm3=0.01,
    )
    metrics = dvh_metrics(samples)

    assert metrics.volume_cm3 == pytest.approx(1.0)
    assert metrics.d98_gy == pytest.approx(np.percentile(samples.dose_gy, 2))
    assert metrics.d95_gy == pytest.approx(np.percentile(samples.dose_gy, 5))
    assert metrics.d50_gy == pytest.approx(50.5)
    assert metrics.d2_gy == pytest.approx(np.percentile(samples.dose_gy, 98))
    assert metrics.mean_gy == pytest.approx(50.5)



class _ScaledPatientSampler:
    def __init__(
        self,
        scale: float,
        uid: str,
        dose_quantity: DoseQuantity = DoseQuantity.DOSE_TO_MEDIUM,
    ) -> None:
        self.scale = scale
        self.frame_of_reference_uid = uid
        self.dose_quantity = dose_quantity

    def sample_patient_points_gy(self, points_patient_mm: np.ndarray) -> np.ndarray:
        points = np.asarray(points_patient_mm, dtype=np.float64)
        # Synthetic reference field in patient coordinates.
        reference = (
            2.0
            + 1e-3 * points[:, 0]
            + 2e-3 * points[:, 1]
            + 3e-3 * points[:, 2]
        )
        return self.scale * reference


def _synthetic_rtdose() -> RtDose:
    uid = "1.2.3"
    frames, rows, columns = 2, 2, 3
    geometry = DoseGeometry(
        rows=rows,
        columns=columns,
        frames=frames,
        pixel_spacing_mm=(2.0, 3.0),
        image_orientation_patient=(1.0, 0.0, 0.0, 0.0, 1.0, 0.0),
        image_position_patient_mm=(10.0, 20.0, 30.0),
        frame_offsets_mm=np.asarray([0.0, 2.5]),
        frame_positions_patient_mm=np.asarray(
            [[10.0, 20.0, 30.0], [10.0, 20.0, 32.5]]
        ),
        frame_offset_mode="RELATIVE",
    )

    dose = np.empty((frames, rows, columns), dtype=np.float64)
    for frame in range(frames):
        for row in range(rows):
            for column in range(columns):
                point = geometry.voxel_center_patient_mm(frame, row, column)
                dose[frame, row, column] = (
                    2.0
                    + 1e-3 * point[0]
                    + 2e-3 * point[1]
                    + 3e-3 * point[2]
                )

    return RtDose(
        dose=dose,
        geometry=geometry,
        dose_units="GY",
        dose_type="PHYSICAL",
        dose_summation_type="PLAN",
        frame_of_reference_uid=uid,
        sop_instance_uid="9.8.7",
        referenced_rtplan_uids=("5.4.3",),
        referenced_beam_numbers=(),
    )


def test_voxelwise_difference_keeps_tps_as_reference_grid() -> None:
    reference = _synthetic_rtdose()
    result = compare_on_rtdose_grid(
        reference,
        _ScaledPatientSampler(1.02, reference.frame_of_reference_uid),
        reference_dose_quantity=DoseQuantity.DOSE_TO_MEDIUM,
        threshold_fraction_of_reference_max=0.0,
    )

    np.testing.assert_allclose(
        result.evaluated_dose_gy,
        1.02 * reference.dose,
        rtol=0.0,
        atol=1e-12,
    )
    np.testing.assert_allclose(
        result.relative_difference_percent,
        2.0,
        rtol=0.0,
        atol=1e-12,
    )
    assert result.summary.mean_relative_difference_percent == pytest.approx(2.0)


def test_voxelwise_difference_rejects_frame_mismatch() -> None:
    reference = _synthetic_rtdose()
    with pytest.raises(ValueError, match="FrameOfReferenceUID"):
        compare_on_rtdose_grid(
            reference,
            _ScaledPatientSampler(1.0, "different"),
            reference_dose_quantity=DoseQuantity.DOSE_TO_MEDIUM,
        )



def test_voxelwise_difference_requires_explicit_tps_dose_quantity() -> None:
    reference = _synthetic_rtdose()

    with pytest.raises(ValueError, match="TPS dose quantity is unknown"):
        compare_on_rtdose_grid(
            reference,
            _ScaledPatientSampler(1.0, reference.frame_of_reference_uid),
            reference_dose_quantity=None,
        )


def test_voxelwise_difference_rejects_dose_quantity_mismatch() -> None:
    reference = _synthetic_rtdose()

    with pytest.raises(ValueError, match="Dose quantity mismatch"):
        compare_on_rtdose_grid(
            reference,
            _ScaledPatientSampler(
                1.0,
                reference.frame_of_reference_uid,
                DoseQuantity.DOSE_TO_WATER,
            ),
            reference_dose_quantity=DoseQuantity.DOSE_TO_MEDIUM,
        )
