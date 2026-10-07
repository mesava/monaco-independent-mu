from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class RoundedLeafTipTangentGeometry:
    """Idealized source-to-rounded-tip tangent geometry.

    This helper models one cylindrical leaf tip in the 2-D leaf-motion plane.
    It is a research geometry primitive, not a production Agility mapper.

    Coordinate convention:
    - source at (x, z) = (0, 0)
    - isocenter plane at z = SAD
    - cylinder centre at (c, CIL)
    - radius = R
    - projected leaf edge x_iso is the source ray tangent to the cylinder.
    """

    sad_cm: float
    radius_cm: float
    cylinder_axis_z_cm: float

    def __post_init__(self) -> None:
        values = np.asarray(
            [self.sad_cm, self.radius_cm, self.cylinder_axis_z_cm],
            dtype=np.float64,
        )
        if not np.all(np.isfinite(values)):
            raise ValueError("Rounded-tip geometry values must be finite.")
        if self.sad_cm <= 0:
            raise ValueError("SAD must be positive.")
        if self.radius_cm <= 0:
            raise ValueError("Leaf-tip radius must be positive.")
        if self.cylinder_axis_z_cm <= self.radius_cm:
            raise ValueError(
                "Cylinder axis must lie farther from the source than its radius."
            )

    def projected_edge_to_cylinder_origin_cm(
        self,
        projected_edge_cm: np.ndarray | float,
        *,
        opening_side: str,
    ) -> np.ndarray:
        """Convert an isocenter-projected edge to cylinder-centre x.

        opening_side describes which side of the circular tip faces the beam
        opening: negative or positive.

        For ray x = m z, m = x_iso / SAD, tangency requires

            abs(c - m*CIL) / sqrt(1 + m^2) = R.

        Hence

            c = m*CIL - R*sqrt(1+m^2)   for negative opening side
            c = m*CIL + R*sqrt(1+m^2)   for positive opening side.
        """

        if opening_side not in {"negative", "positive"}:
            raise ValueError("opening_side must be 'negative' or 'positive'.")

        x_iso = np.asarray(projected_edge_cm, dtype=np.float64)
        if not np.all(np.isfinite(x_iso)):
            raise ValueError("Projected edge coordinates must be finite.")

        m = x_iso / self.sad_cm
        sign = -1.0 if opening_side == "negative" else 1.0
        return (
            m * self.cylinder_axis_z_cm
            + sign * self.radius_cm * np.sqrt(1.0 + m * m)
        )

    def cylinder_origin_to_projected_edge_cm(
        self,
        cylinder_origin_cm: np.ndarray | float,
        *,
        opening_side: str,
    ) -> np.ndarray:
        """Invert the idealized tangent relation.

        The quadratic has two mathematical tangents. The requested opening
        side selects the branch whose signed distance has the physical sign.
        """

        if opening_side not in {"negative", "positive"}:
            raise ValueError("opening_side must be 'negative' or 'positive'.")

        c = np.asarray(cylinder_origin_cm, dtype=np.float64)
        if not np.all(np.isfinite(c)):
            raise ValueError("Cylinder-origin coordinates must be finite.")

        zc = self.cylinder_axis_z_cm
        radius = self.radius_cm
        denominator = zc * zc - radius * radius
        radical = radius * np.sqrt(zc * zc + c * c - radius * radius)

        m1 = (c * zc + radical) / denominator
        m2 = (c * zc - radical) / denominator

        requested_sign = -1.0 if opening_side == "negative" else 1.0
        residual1 = c - m1 * zc
        residual2 = c - m2 * zc

        valid1 = requested_sign * residual1 > 0
        valid2 = requested_sign * residual2 > 0

        if np.any(valid1 == valid2):
            raise ValueError(
                "Cannot select a unique tangent branch for the supplied "
                "cylinder coordinate and opening side."
            )

        slope = np.where(valid1, m1, m2)
        return slope * self.sad_cm

    def tangent_point_cm(
        self,
        projected_edge_cm: np.ndarray | float,
        cylinder_origin_cm: np.ndarray | float,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Return the source-ray tangent point (x, z) on the cylinder."""

        x_iso = np.asarray(projected_edge_cm, dtype=np.float64)
        c = np.asarray(cylinder_origin_cm, dtype=np.float64)
        if not np.all(np.isfinite(x_iso)) or not np.all(np.isfinite(c)):
            raise ValueError("Tangent-point inputs must be finite.")

        m = x_iso / self.sad_cm
        z_tangent = (
            c * m + self.cylinder_axis_z_cm
        ) / (1.0 + m * m)
        x_tangent = m * z_tangent
        return x_tangent, z_tangent

    def tangent_distance_cm(
        self,
        projected_edge_cm: np.ndarray | float,
        cylinder_origin_cm: np.ndarray | float,
    ) -> np.ndarray:
        """Perpendicular distance from cylinder centre to source ray."""

        x_iso = np.asarray(projected_edge_cm, dtype=np.float64)
        c = np.asarray(cylinder_origin_cm, dtype=np.float64)
        m = x_iso / self.sad_cm
        return np.abs(c - m * self.cylinder_axis_z_cm) / np.sqrt(1.0 + m * m)

    def tangent_within_leaf_slab(
        self,
        projected_edge_cm: np.ndarray | float,
        cylinder_origin_cm: np.ndarray | float,
        *,
        zmin_cm: float,
        zmax_cm: float,
        tolerance_cm: float = 1e-9,
    ) -> np.ndarray:
        """Check that the tangent point lies on the physical rounded-tip slab."""

        if not (0.0 < zmin_cm < zmax_cm):
            raise ValueError("Leaf slab requires 0 < zmin_cm < zmax_cm.")
        if tolerance_cm < 0:
            raise ValueError("tolerance_cm must be non-negative.")

        _, z_tangent = self.tangent_point_cm(
            projected_edge_cm,
            cylinder_origin_cm,
        )
        return (
            (z_tangent >= zmin_cm - tolerance_cm)
            & (z_tangent <= zmax_cm + tolerance_cm)
        )


def agility_literature_tangent_geometry(
    *,
    sad_cm: float = 100.0,
) -> RoundedLeafTipTangentGeometry:
    """Research-only Agility values reported in public literature.

    Radius = 17 cm and centre-of-curvature source distance = 34.93 cm.
    These values do not define the complete Agility leaf/body geometry.
    """

    return RoundedLeafTipTangentGeometry(
        sad_cm=sad_cm,
        radius_cm=17.0,
        cylinder_axis_z_cm=34.93,
    )



def leaf_bank_rotation_shift_mm(
    *,
    leaf_thickness_mm: float,
    lbrot_rad: float,
) -> float:
    """Return the magnitude of the field-centre shift caused by LBROT.

    Gholampourkashi et al. used

        shift = leaf_thickness * sin(LBROT) / 2

    to compensate the translation introduced by the Agility leaf-bank
    rotation in the SYNCMLCE model.

    This helper returns only the magnitude.  The sign depends on the chosen
    bank/axis convention and must not be inferred here.
    """

    if not np.isfinite(leaf_thickness_mm) or leaf_thickness_mm <= 0:
        raise ValueError("leaf_thickness_mm must be finite and positive.")
    if not np.isfinite(lbrot_rad):
        raise ValueError("lbrot_rad must be finite.")

    return float(abs(leaf_thickness_mm * np.sin(lbrot_rad) / 2.0))
