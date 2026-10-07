from __future__ import annotations

from dataclasses import dataclass

from .tps_reference import TpsReferenceDataset


@dataclass(frozen=True)
class TpsReferenceComparison:
    energy_model_id: str
    independent_dose_gy: float
    tps_reference_dose_gy: float
    delivered_mu: float
    difference_percent: float


def compare_independent_dose_to_tps_reference(
    *,
    energy_model_id: str,
    independent_dose_gy: float,
    tps_reference: TpsReferenceDataset,
) -> TpsReferenceComparison:
    """Compare an already-normalized independent dose against Monaco.

    This function performs no renormalization and no fitting. The independent
    dose must already be in Gy from an independent absolute calibration.
    """

    if independent_dose_gy <= 0:
        raise ValueError("independent_dose_gy must be positive.")

    reference = tps_reference.energy(energy_model_id)
    difference = 100.0 * (
        independent_dose_gy - reference.dose_gy
    ) / reference.dose_gy

    return TpsReferenceComparison(
        energy_model_id=energy_model_id,
        independent_dose_gy=float(independent_dose_gy),
        tps_reference_dose_gy=reference.dose_gy,
        delivered_mu=tps_reference.geometry.delivered_mu,
        difference_percent=float(difference),
    )
