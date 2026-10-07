import numpy as np
import pytest

from indep_mu.comparison3d import PatientMcDoseGridGy, dvh_metrics
from indep_mu.comparison3d.dvh import StructureDoseSamples
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
