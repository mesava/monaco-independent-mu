import numpy as np
import pytest

from indep_mu.beam_model.iec_coordinates import IsocenterProjection, split_dicom_banks


def test_split_dicom_bank_order_is_preserved() -> None:
    banks = split_dicom_banks(
        [-10.0, -8.0, 8.0, 10.0],
        pair_count=2,
    )

    np.testing.assert_allclose(banks.bank1_mm, [-10.0, -8.0])
    np.testing.assert_allclose(banks.bank2_mm, [8.0, 10.0])


def test_isocenter_projection_round_trip() -> None:
    projection = IsocenterProjection(sad_mm=1000.0)

    physical = projection.to_source_distance_cm(
        np.asarray([-50.0, 50.0]),
        source_distance_cm=40.0,
    )
    np.testing.assert_allclose(physical, [-2.0, 2.0])

    restored = projection.to_isocenter_mm(
        physical,
        source_distance_cm=40.0,
    )
    np.testing.assert_allclose(restored, [-50.0, 50.0])


def test_projection_rejects_invalid_geometry() -> None:
    with pytest.raises(ValueError):
        IsocenterProjection(0.0)
