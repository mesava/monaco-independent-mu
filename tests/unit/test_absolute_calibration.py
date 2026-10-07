import numpy as np
import pytest

from indep_mu.beam_model.commissioning import FieldSize, ReferenceCalibration
from indep_mu.montecarlo.absolute_calibration import AbsoluteMcCalibration
from indep_mu.montecarlo.threedose import ThreeDDose


def _measurement() -> ReferenceCalibration:
    return ReferenceCalibration(
        field=FieldSize(100.0, 100.0),
        ssd_mm=900.0,
        depth_mm=100.0,
        delivered_mu=100.0,
        measured_dose_gy=1.000,
    )


def _mc(value: float) -> ThreeDDose:
    return ThreeDDose(
        x_bound_cm=np.asarray([0.0, 1.0]),
        y_bound_cm=np.asarray([0.0, 1.0]),
        z_bound_cm=np.asarray([0.0, 1.0]),
        dose=np.asarray([[[value]]], dtype=np.float64),
        relative_uncertainty=np.asarray([[[0.01]]], dtype=np.float64),
    )


def test_absolute_calibration_converts_history_to_mu() -> None:
    calibration = AbsoluteMcCalibration(
        energy_model_id="6MV",
        reference_measurement=_measurement(),
        mc_reference_dose_gy_per_history=2.0e-12,
    )

    assert calibration.measured_reference_dose_gy_per_mu == pytest.approx(0.01)
    assert calibration.histories_per_mu == pytest.approx(5.0e9)

    patient = _mc(1.5e-12)
    per_mu = calibration.dose_per_mu(patient)
    assert per_mu.dose[0, 0, 0] == pytest.approx(0.0075)

    beam = calibration.dose_for_beam(patient, beam_mu=200.0)
    assert beam.dose[0, 0, 0] == pytest.approx(1.5)
    assert beam.relative_uncertainty[0, 0, 0] == pytest.approx(0.01)


def test_rejects_nonpositive_beam_mu() -> None:
    calibration = AbsoluteMcCalibration(
        energy_model_id="6MV",
        reference_measurement=_measurement(),
        mc_reference_dose_gy_per_history=2.0e-12,
    )

    with pytest.raises(ValueError):
        calibration.dose_for_beam(_mc(1e-12), beam_mu=0.0)
