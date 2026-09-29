import numpy as np

from indep_mu.patient_model.density import red_to_mass_density


def test_red_to_density_key_points() -> None:
    red = np.asarray([0.0, 0.001, 1.0, 2.335])
    rho = red_to_mass_density(red)

    np.testing.assert_allclose(
        rho,
        [
            0.0,
            (np.sqrt(0.99**2 + 4 * 0.01 * 0.001) - 0.99) / 0.02,
            1.0,
            (2.335 - 0.15) / 0.85,
        ],
        rtol=0.0,
        atol=1e-12,
    )


def test_red_to_density_is_continuous_at_one() -> None:
    eps = 1e-9
    left = red_to_mass_density(np.asarray([1.0 - eps]))[0]
    exact = red_to_mass_density(np.asarray([1.0]))[0]
    right = red_to_mass_density(np.asarray([1.0 + eps]))[0]

    assert abs(left - exact) < 2e-9
    assert abs(right - exact) < 2e-9
