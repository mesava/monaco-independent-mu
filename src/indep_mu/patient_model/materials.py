from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable


@dataclass(frozen=True)
class MaterialThreshold:
    material_name: str
    mass_density_g_cm3: float


# Raw patient LUT transcribed from the supplied Monaco configuration.
# Interpretation/interpolation between different adjacent materials is handled
# separately; this table is kept lossless and ordered.
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
