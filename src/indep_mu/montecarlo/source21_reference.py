from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class CoplanarHfsReferenceAngles:
    """Published coplanar HFS reference transform, not a production mapper."""

    theta_deg: float
    phi_deg: float
    phicol_deg: float


def coplanar_hfs_reference_angles(
    *,
    gantry_deg: float,
    collimator_deg: float,
    patient_support_deg: float,
    couch_tolerance_deg: float = 1e-6,
) -> CoplanarHfsReferenceAngles:
    """Return the simple published DICOM->DOSXYZnrc coplanar reference.

    This implements the commonly cited coplanar HFS relation:

        theta   = 90 deg
        phi     = -90 deg + gantry
        phicol  = -90 deg - collimator

    It is intentionally *not* an implementation of Source21OrientationMapper.
    The function exists only to provide regression vectors while the general
    non-coplanar rotation-matrix transform is derived and validated.

    A non-zero patient-support/couch angle is rejected.
    """

    couch = float(patient_support_deg) % 360.0
    couch_distance = min(abs(couch), abs(couch - 360.0))
    if couch_distance > couch_tolerance_deg:
        raise ValueError(
            "Coplanar HFS reference transform is valid only for couch angle 0."
        )

    return CoplanarHfsReferenceAngles(
        theta_deg=90.0,
        phi_deg=(-90.0 + float(gantry_deg)) % 360.0,
        phicol_deg=(-90.0 - float(collimator_deg)) % 360.0,
    )
