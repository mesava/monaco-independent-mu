from __future__ import annotations

import numpy as np


def red_to_mass_density(red_values: np.ndarray) -> np.ndarray:
    """Convert relative electron density (RED) to mass density in g/cm3.

    Implements the piecewise Monaco patient-model relation documented by
    Elekta. The function is vectorized and continuous at RED = 1.
    """
    red = np.asarray(red_values, dtype=np.float64)
    if np.any(~np.isfinite(red)):
        raise ValueError("RED values must be finite.")

    rho = np.empty_like(red, dtype=np.float64)

    non_positive = red <= 0.0
    below_water = (red > 0.0) & (red < 1.0)
    at_or_above_water = red >= 1.0

    rho[non_positive] = 0.0
    rho[below_water] = (
        np.sqrt(0.99**2 + 4.0 * 0.01 * red[below_water]) - 0.99
    ) / (2.0 * 0.01)
    rho[at_or_above_water] = (red[at_or_above_water] - 0.15) / 0.85

    return rho
