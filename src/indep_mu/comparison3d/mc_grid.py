from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from indep_mu.montecarlo.threedose import ThreeDDose
from indep_mu.mu_check.dosegrid import RectilinearDoseGrid, trilinear_sample
from indep_mu.patient_model.egsphant import EgsphantGeometry

from .dose_quantity import DoseQuantity


@dataclass(frozen=True)
class PatientMcDoseGridGy:
    """Absolute MC dose grid tied to the DICOM patient coordinate frame.

    The dose must already be independently normalized to Gy.
    """

    dose: ThreeDDose
    phantom_geometry: EgsphantGeometry
    frame_of_reference_uid: str
    dose_quantity: DoseQuantity

    def __post_init__(self) -> None:
        if not self.frame_of_reference_uid:
            raise ValueError("frame_of_reference_uid is required.")

        np.testing.assert_allclose(
            np.asarray(self.dose.x_bound_cm, dtype=np.float64),
            np.asarray(self.phantom_geometry.x_bound_cm, dtype=np.float64),
            rtol=0.0,
            atol=1e-8,
            err_msg="MC x boundaries do not match egsphant geometry.",
        )
        np.testing.assert_allclose(
            np.asarray(self.dose.y_bound_cm, dtype=np.float64),
            np.asarray(self.phantom_geometry.y_bound_cm, dtype=np.float64),
            rtol=0.0,
            atol=1e-8,
            err_msg="MC y boundaries do not match egsphant geometry.",
        )
        np.testing.assert_allclose(
            np.asarray(self.dose.z_bound_cm, dtype=np.float64),
            np.asarray(self.phantom_geometry.z_bound_cm, dtype=np.float64),
            rtol=0.0,
            atol=1e-8,
            err_msg="MC z boundaries do not match egsphant geometry.",
        )

    @property
    def rectilinear_grid(self) -> RectilinearDoseGrid:
        x = 0.5 * (self.dose.x_bound_cm[:-1] + self.dose.x_bound_cm[1:])
        y = 0.5 * (self.dose.y_bound_cm[:-1] + self.dose.y_bound_cm[1:])
        z = 0.5 * (self.dose.z_bound_cm[:-1] + self.dose.z_bound_cm[1:])
        return RectilinearDoseGrid(
            x_cm=np.asarray(x, dtype=np.float64),
            y_cm=np.asarray(y, dtype=np.float64),
            z_cm=np.asarray(z, dtype=np.float64),
            dose=np.asarray(self.dose.dose, dtype=np.float64),
        )

    def patient_to_mc_cm(self, points_patient_mm: np.ndarray) -> np.ndarray:
        """Project DICOM patient coordinates onto the egsphant axis basis."""

        points = np.asarray(points_patient_mm, dtype=np.float64)
        if points.ndim == 1:
            points = points.reshape(1, -1)
        if points.ndim != 2 or points.shape[1] != 3:
            raise ValueError("Patient points must have shape (N, 3).")
        if not np.all(np.isfinite(points)):
            raise ValueError("Patient points must be finite.")

        axes = np.vstack(
            [
                self.phantom_geometry.x_direction_patient,
                self.phantom_geometry.y_direction_patient,
                self.phantom_geometry.z_direction_patient,
            ]
        )
        return (points @ axes.T) / 10.0

    def sample_patient_points_gy(self, points_patient_mm: np.ndarray) -> np.ndarray:
        local = self.patient_to_mc_cm(points_patient_mm)
        return trilinear_sample(self.rectilinear_grid, local)
