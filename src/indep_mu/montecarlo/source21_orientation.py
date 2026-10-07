from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .source21 import Source21Angles


def _angular_distance_deg(a: float, b: float) -> float:
    delta = (float(a) - float(b)) % 360.0
    return min(delta, 360.0 - delta)


@dataclass(frozen=True)
class ZhanHfsSource21OrientationMapper:
    """DICOM IEC angles -> DOSXYZnrc source-21 angles for HFS geometry.

    The implementation follows the compact transform published by Zhan,
    Jiang and Osei (Phys Med Biol 2012) and the reference implementation
    released by Lixin Zhan.

    Exact singular combinations couch=90/270 and gantry=90/270 are rejected
    instead of silently perturbing the clinical angles.  At those points the
    spherical azimuth/collimator decomposition is singular even though the
    physical beam orientation remains well-defined.
    """

    patient_position: str = "HFS"
    singularity_tolerance_deg: float = 1e-8

    def __post_init__(self) -> None:
        if self.patient_position.upper() != "HFS":
            raise ValueError(
                "ZhanHfsSource21OrientationMapper is validated only for HFS."
            )
        if self.singularity_tolerance_deg <= 0:
            raise ValueError("singularity_tolerance_deg must be positive.")

    def _is_singular(self, gantry_deg: float, couch_deg: float) -> bool:
        gantry_cardinal = any(
            _angular_distance_deg(gantry_deg, target)
            <= self.singularity_tolerance_deg
            for target in (90.0, 270.0)
        )
        couch_cardinal = any(
            _angular_distance_deg(couch_deg, target)
            <= self.singularity_tolerance_deg
            for target in (90.0, 270.0)
        )
        return gantry_cardinal and couch_cardinal

    def map_angles(
        self,
        *,
        gantry_deg: float,
        collimator_deg: float,
        patient_support_deg: float,
    ) -> Source21Angles:
        values = np.asarray(
            [gantry_deg, collimator_deg, patient_support_deg],
            dtype=np.float64,
        )
        if not np.all(np.isfinite(values)):
            raise ValueError("DICOM beam angles must be finite.")

        if self._is_singular(gantry_deg, patient_support_deg):
            raise ValueError(
                "DICOM->DOSXYZnrc transform is singular for simultaneous "
                "gantry and couch angles at 90/270 degrees.  The project "
                "rejects this control point instead of perturbing angles."
            )

        gamma = np.deg2rad(float(gantry_deg))
        rho = np.deg2rad(float(patient_support_deg))
        col = np.deg2rad(float(collimator_deg))

        sin_gamma = np.sin(gamma)
        cos_gamma = np.cos(gamma)
        sin_rho = np.sin(rho)
        cos_rho = np.cos(rho)

        sgsr = sin_gamma * sin_rho
        sgcr = sin_gamma * cos_rho

        # Clamp only round-off excursions outside [-1, 1].
        theta = np.arccos(np.clip(-sgsr, -1.0, 1.0))
        phi = np.arctan2(-cos_gamma, sgcr)

        couch_angle_in_collimator_plane = np.arctan2(
            -sin_rho * cos_gamma,
            cos_rho,
        )

        # Zhan reference implementation:
        # phicol = (col - pi/2) + couch_projection
        # followed by the BEAMnrc phase-space -> DOSXYZnrc transform.
        phicol = (
            np.pi
            - (
                (col - np.pi / 2.0)
                + couch_angle_in_collimator_plane
            )
        )

        return Source21Angles(
            theta_deg=float(np.rad2deg(theta)),
            phi_deg=float(np.mod(np.rad2deg(phi), 360.0)),
            phicol_deg=float(np.mod(np.rad2deg(phicol), 360.0)),
        )
