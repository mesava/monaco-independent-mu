import numpy as np
import pytest

from indep_mu.patient_model.hu_red import DRT120KV


def test_calibration_points_are_exact() -> None:
    calculated = DRT120KV.to_red(DRT120KV.hu)
    np.testing.assert_allclose(calculated, DRT120KV.red, rtol=0.0, atol=1e-12)


def test_linear_interpolation() -> None:
    hu = np.asarray([(-1000 + -719) / 2])
    expected = np.asarray([(0.001 + 0.227) / 2])
    np.testing.assert_allclose(DRT120KV.to_red(hu), expected)


def test_out_of_range_is_rejected_by_default() -> None:
    with pytest.raises(ValueError):
        DRT120KV.to_red(np.asarray([-1001.0]))


def test_out_of_range_can_be_clipped_explicitly() -> None:
    values = DRT120KV.to_red(np.asarray([-2000.0, 3000.0]), out_of_range="clip")
    np.testing.assert_allclose(values, [0.001, 2.335])
