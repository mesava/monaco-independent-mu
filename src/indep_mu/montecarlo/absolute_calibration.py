from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from indep_mu.beam_model.commissioning import ReferenceCalibration

from .threedose import ThreeDDose


@dataclass(frozen=True)
class AbsoluteMcCalibration:
    """Independent conversion from MC dose/history to clinical Gy/MU.

    The reference MC calculation and the patient calculation must use the same
    source normalization convention for this factor to be transferable.
    """

    energy_model_id: str
    reference_measurement: ReferenceCalibration
    mc_reference_dose_gy_per_history: float

    def __post_init__(self) -> None:
        if not self.energy_model_id:
            raise ValueError("energy_model_id is required.")
        if (
            not np.isfinite(self.mc_reference_dose_gy_per_history)
            or self.mc_reference_dose_gy_per_history <= 0
        ):
            raise ValueError(
                "mc_reference_dose_gy_per_history must be finite and positive."
            )

    @property
    def measured_reference_dose_gy_per_mu(self) -> float:
        return self.reference_measurement.dose_gy_per_mu

    @property
    def histories_per_mu(self) -> float:
        """Equivalent source histories per MU for the commissioned MC model."""

        return (
            self.measured_reference_dose_gy_per_mu
            / self.mc_reference_dose_gy_per_history
        )

    def dose_per_mu(self, mc_dose_per_history: ThreeDDose) -> ThreeDDose:
        """Convert an MC field from Gy/history to Gy/MU."""

        return mc_dose_per_history.scaled(self.histories_per_mu)

    def dose_for_beam(
        self,
        mc_dose_per_history: ThreeDDose,
        *,
        beam_mu: float,
    ) -> ThreeDDose:
        """Convert an MC field from Gy/history to absolute beam dose in Gy."""

        if not np.isfinite(beam_mu) or beam_mu <= 0:
            raise ValueError("beam_mu must be finite and positive.")
        return mc_dose_per_history.scaled(self.histories_per_mu * beam_mu)


def assert_reference_geometry_matches(
    a: ReferenceCalibration,
    b: ReferenceCalibration,
    *,
    tolerance_mm: float = 1e-6,
) -> None:
    """Require two reference calibrations to describe the same geometry."""

    if abs(a.field.x_mm - b.field.x_mm) > tolerance_mm:
        raise ValueError("Reference field x dimension does not match.")
    if abs(a.field.y_mm - b.field.y_mm) > tolerance_mm:
        raise ValueError("Reference field y dimension does not match.")
    if abs(a.ssd_mm - b.ssd_mm) > tolerance_mm:
        raise ValueError("Reference SSD does not match.")
    if abs(a.depth_mm - b.depth_mm) > tolerance_mm:
        raise ValueError("Reference depth does not match.")
