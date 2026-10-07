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


@dataclass(frozen=True)
class MaterialBlendArrays:
    """Compact voxel-wise representation of Monaco Patient-LUT interpolation."""

    material_names: tuple[str, ...]
    lower_material_id: np.ndarray
    upper_material_id: np.ndarray
    upper_fraction: np.ndarray


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


def _material_names_and_entry_ids(
    entries: tuple[MaterialThreshold, ...],
) -> tuple[tuple[str, ...], np.ndarray]:
    names: list[str] = []
    ids: list[int] = []
    for entry in entries:
        if entry.material_name not in names:
            names.append(entry.material_name)
        ids.append(names.index(entry.material_name))
    return tuple(names), np.asarray(ids, dtype=np.uint8)


def material_blend_arrays_for_density(
    mass_density_g_cm3: np.ndarray,
    entries: tuple[MaterialThreshold, ...] = PATIENT_LUT,
) -> MaterialBlendArrays:
    """Vectorized Monaco-style material interpolation.

    Equal density thresholds represent an abrupt transition. Exact ties use the
    right-most threshold, which is important at 3, 5 and 9 g/cm3 in the supplied
    Patient LUT.
    """

    validate_lut(entries)
    rho = np.asarray(mass_density_g_cm3, dtype=np.float64)
    if np.any(~np.isfinite(rho)):
        raise ValueError("Mass density values must be finite.")

    thresholds = np.asarray(
        [entry.mass_density_g_cm3 for entry in entries], dtype=np.float64
    )
    if np.any(rho < thresholds[0]) or np.any(rho > thresholds[-1]):
        raise ValueError(
            f"Mass density outside LUT range [{thresholds[0]}, {thresholds[-1]}] g/cm3."
        )

    material_names, entry_ids = _material_names_and_entry_ids(entries)

    upper_index = np.searchsorted(thresholds, rho, side="right")
    at_last = upper_index >= len(entries)

    upper_safe = np.minimum(upper_index, len(entries) - 1)
    lower_safe = np.maximum(upper_safe - 1, 0)

    lower_id = entry_ids[lower_safe].astype(np.uint8, copy=False)
    upper_id = entry_ids[upper_safe].astype(np.uint8, copy=False)

    lower_density = thresholds[lower_safe]
    upper_density = thresholds[upper_safe]
    denominator = upper_density - lower_density

    fraction = np.zeros(rho.shape, dtype=np.float64)
    interpolated = (~at_last) & (denominator > 0.0)
    fraction[interpolated] = (
        (rho[interpolated] - lower_density[interpolated])
        / denominator[interpolated]
    )

    # Repeated names represent one pure composition over a density interval.
    same_material = lower_id == upper_id
    fraction[same_material | at_last] = 0.0

    return MaterialBlendArrays(
        material_names=material_names,
        lower_material_id=lower_id,
        upper_material_id=upper_id,
        upper_fraction=fraction.astype(np.float32),
    )


def material_blend_for_density(
    mass_density_g_cm3: float,
    entries: tuple[MaterialThreshold, ...] = PATIENT_LUT,
) -> MaterialBlend:
    """Return Monaco-style material interpolation for one mass density."""

    field = material_blend_arrays_for_density(
        np.asarray([mass_density_g_cm3], dtype=np.float64),
        entries,
    )
    lower_id = int(field.lower_material_id[0])
    upper_id = int(field.upper_material_id[0])

    thresholds = np.asarray(
        [entry.mass_density_g_cm3 for entry in entries], dtype=np.float64
    )
    rho = float(mass_density_g_cm3)
    upper_index = int(np.searchsorted(thresholds, rho, side="right"))

    if upper_index >= len(entries):
        lower_density = upper_density = float(thresholds[-1])
    else:
        lower_index = max(0, upper_index - 1)
        lower_density = float(thresholds[lower_index])
        upper_density = float(thresholds[upper_index])

    return MaterialBlend(
        lower_material=field.material_names[lower_id],
        upper_material=field.material_names[upper_id],
        lower_density_g_cm3=lower_density,
        upper_density_g_cm3=upper_density,
        upper_fraction=float(field.upper_fraction[0]),
    )
