import numpy as np
import pytest

from indep_mu.mu_check import (
    RectilinearDoseGrid,
    SphericalRoi,
    check_beam_equivalent_mu,
    spherical_roi_mean,
    trilinear_sample,
)


def _linear_grid() -> RectilinearDoseGrid:
    x = np.linspace(-1.0, 1.0, 9)
    y = np.linspace(-1.0, 1.0, 9)
    z = np.linspace(-1.0, 1.0, 9)
    zz, yy, xx = np.meshgrid(z, y, x, indexing="ij")
    dose = 2.0 + 0.3 * xx - 0.2 * yy + 0.1 * zz
    return RectilinearDoseGrid(x_cm=x, y_cm=y, z_cm=z, dose=dose)


def test_trilinear_interpolation_is_exact_for_linear_field() -> None:
    grid = _linear_grid()
    point = np.asarray([[0.13, -0.27, 0.41]])
    expected = 2.0 + 0.3 * 0.13 - 0.2 * (-0.27) + 0.1 * 0.41
    assert trilinear_sample(grid, point)[0] == pytest.approx(expected)


def test_symmetric_spherical_mean_of_linear_field_equals_centre() -> None:
    grid = _linear_grid()
    roi = SphericalRoi("test", (0.2, -0.1, 0.3), radius_cm=0.25)
    result = spherical_roi_mean(grid, roi, samples_per_axis=11)

    expected = 2.0 + 0.3 * 0.2 - 0.2 * (-0.1) + 0.1 * 0.3
    assert result.mean_dose == pytest.approx(expected, abs=1e-12)
    assert result.sample_count > 100


def test_equivalent_mu_multi_roi_check() -> None:
    result = check_beam_equivalent_mu(
        beam_number=1,
        tps_mu=200.0,
        roi_doses_gy=(
            ("ROI1", 1.00, 0.98),
            ("ROI2", 2.00, 2.02),
            ("ROI3", 1.50, 1.50),
        ),
        tolerance_percent=5.0,
    )

    assert result.rois[0].equivalent_mu == pytest.approx(204.081632653)
    assert result.rois[0].difference_percent == pytest.approx(2.0408163265)
    assert result.rois[1].difference_percent == pytest.approx(-0.9900990099)
    assert result.within_tolerance
    assert result.roi_spread_percent_of_tps_mu > 0


def test_equivalent_mu_flags_tolerance_exceedance() -> None:
    result = check_beam_equivalent_mu(
        beam_number=2,
        tps_mu=100.0,
        roi_doses_gy=(("ROI1", 1.0, 0.90),),
        tolerance_percent=5.0,
    )
    assert result.max_abs_difference_percent == pytest.approx(11.111111111)
    assert not result.within_tolerance
