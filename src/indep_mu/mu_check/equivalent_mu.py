from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class RoiEquivalentMu:
    roi_name: str
    tps_dose_gy: float
    mc_dose_gy_at_tps_mu: float
    equivalent_mu: float
    difference_percent: float


@dataclass(frozen=True)
class BeamEquivalentMuCheck:
    beam_number: int
    tps_mu: float
    rois: tuple[RoiEquivalentMu, ...]
    tolerance_percent: float

    @property
    def max_abs_difference_percent(self) -> float:
        return max(abs(item.difference_percent) for item in self.rois)

    @property
    def roi_spread_percent_of_tps_mu(self) -> float:
        values = [item.equivalent_mu for item in self.rois]
        return 100.0 * (max(values) - min(values)) / self.tps_mu

    @property
    def within_tolerance(self) -> bool:
        return self.max_abs_difference_percent <= self.tolerance_percent


def equivalent_mu_from_delivered_doses(
    *,
    roi_name: str,
    tps_mu: float,
    tps_dose_gy: float,
    mc_dose_gy_at_tps_mu: float,
) -> RoiEquivalentMu:
    """Equivalent MU consistency check for a TPS-delivered beam.

    MU_equiv = MU_TPS * D_TPS / D_MC(MU_TPS)

    This is intentionally named *equivalent MU*. For an optimized VMAT beam it
    does not independently reconstruct the optimizer's MU from prescription.
    It tests the independent MC dose-per-MU against the TPS beam dose.
    """

    values = (tps_mu, tps_dose_gy, mc_dose_gy_at_tps_mu)
    if any(not np.isfinite(value) or value <= 0 for value in values):
        raise ValueError("MU and ROI doses must be finite and positive.")

    equivalent = tps_mu * tps_dose_gy / mc_dose_gy_at_tps_mu
    difference = 100.0 * (equivalent - tps_mu) / tps_mu

    return RoiEquivalentMu(
        roi_name=roi_name,
        tps_dose_gy=float(tps_dose_gy),
        mc_dose_gy_at_tps_mu=float(mc_dose_gy_at_tps_mu),
        equivalent_mu=float(equivalent),
        difference_percent=float(difference),
    )


def check_beam_equivalent_mu(
    *,
    beam_number: int,
    tps_mu: float,
    roi_doses_gy: tuple[tuple[str, float, float], ...],
    tolerance_percent: float,
) -> BeamEquivalentMuCheck:
    """Build a multi-ROI beam check.

    roi_doses_gy entries are (roi_name, TPS beam dose, MC beam dose at TPS MU).
    """

    if not roi_doses_gy:
        raise ValueError("At least one ROI is required.")
    if not np.isfinite(tolerance_percent) or tolerance_percent <= 0:
        raise ValueError("tolerance_percent must be finite and positive.")

    rois = tuple(
        equivalent_mu_from_delivered_doses(
            roi_name=name,
            tps_mu=tps_mu,
            tps_dose_gy=tps_dose,
            mc_dose_gy_at_tps_mu=mc_dose,
        )
        for name, tps_dose, mc_dose in roi_doses_gy
    )

    return BeamEquivalentMuCheck(
        beam_number=int(beam_number),
        tps_mu=float(tps_mu),
        rois=rois,
        tolerance_percent=float(tolerance_percent),
    )
