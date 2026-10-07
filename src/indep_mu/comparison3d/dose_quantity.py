from __future__ import annotations

from enum import Enum


class DoseQuantity(str, Enum):
    """Physical dose quantity used for TPS/MC comparison."""

    DOSE_TO_MEDIUM = "DOSE_TO_MEDIUM"
    DOSE_TO_WATER = "DOSE_TO_WATER"


def require_matching_dose_quantity(
    reference: DoseQuantity | None,
    evaluated: DoseQuantity | None,
) -> DoseQuantity:
    """Require an explicit, identical physical dose quantity.

    Standard RTDOSE tags identify units/type/summation but do not prove whether
    a Monaco Monte Carlo plan was reported as dose-to-medium or converted to
    dose-to-water.  Therefore the TPS quantity must be supplied from verified
    calculation metadata/workflow rather than inferred from DoseUnits=GY.
    """

    if reference is None:
        raise ValueError(
            "TPS dose quantity is unknown. Explicitly declare verified "
            "DOSE_TO_MEDIUM or DOSE_TO_WATER before 3D comparison."
        )
    if evaluated is None:
        raise ValueError(
            "Independent MC dose quantity is unknown. Explicitly declare the "
            "scored/converted quantity before 3D comparison."
        )
    if reference != evaluated:
        raise ValueError(
            "Dose quantity mismatch: "
            f"TPS={reference.value}, independent={evaluated.value}."
        )
    return reference
