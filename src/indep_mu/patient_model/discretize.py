from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .materials import MaterialBlendArrays

# Standard ctcreate .egsphant encoding contains 62 characters and uses
# encoding(medium_index + 1), leaving 61 usable 1-based medium numbers.
EGSPHANT_MAX_MEDIA = 61


@dataclass(frozen=True)
class DiscreteMedium:
    egsphant_name: str
    material_a: str
    material_b: str
    fraction_b: float

    @property
    def is_pure(self) -> bool:
        return self.material_a == self.material_b or self.fraction_b == 0.0


@dataclass(frozen=True)
class DiscreteMaterialField:
    """Finite-medium representation suitable for legacy .egsphant."""

    medium_index: np.ndarray  # 1-based medium number
    media: tuple[DiscreteMedium, ...]
    mixture_bins: int


def _medium_name(a: int, b: int, bin_index: int) -> str:
    if a == b:
        return f"MAT{a:02d}"
    return f"MIX{a:02d}{b:02d}{bin_index:02d}"


def discretize_material_blends(
    blends: MaterialBlendArrays,
    *,
    mixture_bins: int,
    max_media: int = EGSPHANT_MAX_MEDIA,
) -> DiscreteMaterialField:
    """Quantize continuous Monaco material mixtures into finite EGSnrc media.

    mixture_bins is the number of equal intervals between the two pure
    endpoint compositions. For example, mixture_bins=19 gives 18 interior
    compositions per unique material pair.

    No clinical default is chosen deliberately: this parameter must be fixed
    by sensitivity/validation testing.
    """

    if mixture_bins < 1:
        raise ValueError("mixture_bins must be >= 1.")
    if max_media < 1:
        raise ValueError("max_media must be >= 1.")

    lower = np.asarray(blends.lower_material_id, dtype=np.int32)
    upper = np.asarray(blends.upper_material_id, dtype=np.int32)
    fraction = np.asarray(blends.upper_fraction, dtype=np.float64)

    if lower.shape != upper.shape or lower.shape != fraction.shape:
        raise ValueError("Material blend arrays must have identical shapes.")
    if np.any((fraction < 0.0) | (fraction > 1.0)):
        raise ValueError("Material mixture fractions must be in [0, 1].")

    # Canonicalize pair orientation so A<->B and B<->A produce the same
    # physical discrete medium set.
    a = np.minimum(lower, upper)
    b = np.maximum(lower, upper)
    fraction_b = np.where(lower <= upper, fraction, 1.0 - fraction)

    bins = np.rint(fraction_b * mixture_bins).astype(np.int32)
    bins = np.clip(bins, 0, mixture_bins)

    # Collapse endpoints to pure materials.
    pure_a = bins == 0
    pure_b = bins == mixture_bins
    same = a == b

    b = b.copy()
    bins = bins.copy()

    b[pure_a | same] = a[pure_a | same]
    bins[pure_a | same] = 0

    a[pure_b & ~same] = b[pure_b & ~same]
    b[pure_b & ~same] = a[pure_b & ~same]
    bins[pure_b & ~same] = 0

    material_count = len(blends.material_names)
    key_base = material_count
    key = ((a * key_base + b) * (mixture_bins + 1) + bins).astype(np.int64)

    unique_keys, inverse = np.unique(key, return_inverse=True)

    media: list[DiscreteMedium] = []
    for encoded in unique_keys:
        bin_index = int(encoded % (mixture_bins + 1))
        pair_code = int(encoded // (mixture_bins + 1))
        b_id = pair_code % key_base
        a_id = pair_code // key_base

        if a_id == b_id:
            fraction_value = 0.0
        else:
            fraction_value = bin_index / mixture_bins

        media.append(
            DiscreteMedium(
                egsphant_name=_medium_name(a_id, b_id, bin_index),
                material_a=blends.material_names[a_id],
                material_b=blends.material_names[b_id],
                fraction_b=float(fraction_value),
            )
        )

    if len(media) > max_media:
        raise ValueError(
            f"Discretization requires {len(media)} media, exceeding "
            f"egsphant limit {max_media}. Reduce mixture_bins or material set."
        )

    medium_index = (inverse.reshape(lower.shape) + 1).astype(np.uint8)

    return DiscreteMaterialField(
        medium_index=medium_index,
        media=tuple(media),
        mixture_bins=mixture_bins,
    )
