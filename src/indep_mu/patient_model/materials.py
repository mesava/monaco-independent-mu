from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import numpy as np


@dataclass(frozen=True)
class MaterialThreshold:
    material_name: str
    mass_density_g_cm3: float


@dataclass(frozen=True)
class MaterialBlend:
    lower_material: str
    upper_material: str
    lower_density_g_cm3: float
    upper_density_g_cm3: float
    upper_fraction: float

    @property
    def is_pure(self) -> bool:
        return self.lower_material == self.upper_material or self.upper_fraction == 0.0


# Ordered transcription of the supplied Monaco Patient LUT.
PATIENT_LUT: tuple[MaterialThreshold, ...] = (
    MaterialThreshold("DryAir", 0.0),
    MaterialThreshold("DryAir", 0.002),
    MaterialThreshold("MuscleSkeletalIcrp", 0.82),
    MaterialThreshold("AdiposeTissueICRP", 0.91),
    MaterialThreshold("AdiposeTissueICRP", 0.95),
    MaterialThreshold("MuscleSkeletalIcrp", 1.04),
    MaterialThreshold("MuscleSkeletalIcrp", 1.08),
    MaterialThreshold("BoneCorticalIcrp", 1.85),
    MaterialThreshold("BoneCorticalIcrp", 3.0),
    MaterialThreshold("Titanium", 3.0),
    MaterialThreshold("Titanium", 5.0),
    MaterialThreshold("StainlessSteel316", 5.0),
    MaterialThreshold("StainlessSteel316", 9.0),
    MaterialThreshold("Lead", 9.0),
    MaterialThreshold("Lead", 12.0),
)


def validate_lut(entries: Iterable[MaterialThreshold]) -> None:
    items = tuple(entries)
    if not items:
        raise ValueError("Material LUT must not be empty.")
    densities = [entry.mass_density_g_cm3 for entry in items]
    if any(b < a for a, b in zip(densities, densities[1:])):
        raise ValueError("Material LUT density thresholds must be non-decreasing.")


def material_blend_for_density(
    mass_density_g_cm3: float,
    entries: tuple[MaterialThreshold, ...] = PATIENT_LUT,
) -> MaterialBlend:
    """Return Monaco-style material interpolation for one mass density.

    Patient and Phantom LUTs interpolate between the two adjacent LUT
    materials. Repeated identical material names therefore represent a pure
    material over that density interval. Equal density thresholds encode an
    abrupt material transition; exact ties use the right-most threshold.
    """
    validate_lut(entries)
    rho = float(mass_density_g_cm3)
    if not np.isfinite(rho):
        raise ValueError("Mass density must be finite.")

    thresholds = np.asarray(
        [entry.mass_density_g_cm3 for entry in entries], dtype=np.float64
    )
    if rho < thresholds[0] or rho > thresholds[-1]:
        raise ValueError(
            f"Mass density {rho} g/cm3 outside LUT range "
            f"[{thresholds[0]}, {thresholds[-1]}]."
        )

    upper_index = int(np.searchsorted(thresholds, rho, side="right"))
    if upper_index >= len(entries):
        item = entries[-1]
        return MaterialBlend(
            item.material_name, item.material_name,
            item.mass_density_g_cm3, item.mass_density_g_cm3, 0.0
        )

    lower_index = max(0, upper_index - 1)
    lower = entries[lower_index]
    upper = entries[upper_index]

    if upper.mass_density_g_cm3 == lower.mass_density_g_cm3:
        fraction = 0.0
    else:
        fraction = (
            (rho - lower.mass_density_g_cm3)
            / (upper.mass_density_g_cm3 - lower.mass_density_g_cm3)
        )

    return MaterialBlend(
        lower_material=lower.material_name,
        upper_material=upper.material_name,
        lower_density_g_cm3=lower.mass_density_g_cm3,
        upper_density_g_cm3=upper.mass_density_g_cm3,
        upper_fraction=float(fraction),
    )
