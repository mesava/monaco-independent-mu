from pathlib import Path

import numpy as np
import pytest

from indep_mu.montecarlo.threedose import ThreeDDose, read_3ddose, write_3ddose


def _dose() -> ThreeDDose:
    dose = np.arange(12, dtype=np.float64).reshape((2, 2, 3)) / 10.0
    uncertainty = np.full((2, 2, 3), 0.02)
    return ThreeDDose(
        x_bound_cm=np.asarray([-1.5, -0.5, 0.5, 1.5]),
        y_bound_cm=np.asarray([-1.0, 0.0, 1.0]),
        z_bound_cm=np.asarray([0.0, 1.0, 2.0]),
        dose=dose,
        relative_uncertainty=uncertainty,
    )


def test_3ddose_round_trip_preserves_xyz_order(tmp_path: Path) -> None:
    original = _dose()
    path = tmp_path / "test.3ddose"

    write_3ddose(path, original)
    restored = read_3ddose(path)

    assert restored.shape == (2, 2, 3)
    np.testing.assert_allclose(restored.x_bound_cm, original.x_bound_cm)
    np.testing.assert_allclose(restored.y_bound_cm, original.y_bound_cm)
    np.testing.assert_allclose(restored.z_bound_cm, original.z_bound_cm)
    np.testing.assert_allclose(restored.dose, original.dose)
    np.testing.assert_allclose(
        restored.relative_uncertainty,
        original.relative_uncertainty,
    )


def test_scaling_does_not_change_relative_uncertainty() -> None:
    original = _dose()
    scaled = original.scaled(2.5)

    np.testing.assert_allclose(scaled.dose, original.dose * 2.5)
    np.testing.assert_allclose(
        scaled.relative_uncertainty,
        original.relative_uncertainty,
    )


def test_rejects_invalid_boundary_shape() -> None:
    with pytest.raises(ValueError, match="Dose shape"):
        ThreeDDose(
            x_bound_cm=np.asarray([0.0, 1.0]),
            y_bound_cm=np.asarray([0.0, 1.0]),
            z_bound_cm=np.asarray([0.0, 1.0]),
            dose=np.zeros((1, 1, 2)),
            relative_uncertainty=np.zeros((1, 1, 2)),
        )
