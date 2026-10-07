from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .compositions import BASE_COMPOSITIONS
from .discretize import DiscreteMaterialField, discretize_material_blends
from .materials import MaterialBlendArrays
from .pegsless import build_pegsless_media_set


@dataclass(frozen=True)
class DiscretizationSensitivityResult:
    mixture_bins: int
    feasible: bool
    media_count: int | None
    voxel_count: int
    mean_abs_fraction_error: float | None
    p95_abs_fraction_error: float | None
    max_abs_fraction_error: float | None
    rhof_min: float | None
    rhof_max: float | None
    rhof_p95_abs_percent: float | None
    rhof_over_20_percent_fraction: float | None
    error: str | None = None


def _canonical_blend(
    blends: MaterialBlendArrays,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    lower = np.asarray(blends.lower_material_id, dtype=np.int32)
    upper = np.asarray(blends.upper_material_id, dtype=np.int32)
    fraction = np.asarray(blends.upper_fraction, dtype=np.float64)

    a = np.minimum(lower, upper)
    b = np.maximum(lower, upper)
    fraction_b = np.where(lower <= upper, fraction, 1.0 - fraction)
    fraction_b = np.where(a == b, 0.0, fraction_b)
    return a, b, fraction_b


def _quantized_fraction_field(
    blends: MaterialBlendArrays,
    discrete: DiscreteMaterialField,
) -> np.ndarray:
    a, b, _ = _canonical_blend(blends)
    names = blends.material_names
    name_to_id = {name: index for index, name in enumerate(names)}

    quantized = np.empty(a.shape, dtype=np.float64)
    medium_index = np.asarray(discrete.medium_index, dtype=np.int32)

    for index, medium in enumerate(discrete.media, start=1):
        voxels = medium_index == index
        if not np.any(voxels):
            continue

        med_a = name_to_id[medium.material_a]
        med_b = name_to_id[medium.material_b]

        local_a = a[voxels]
        local_b = b[voxels]

        if med_a == med_b:
            is_original_a = med_a == local_a
            is_original_b = med_a == local_b
            if np.any(~(is_original_a | is_original_b)):
                raise ValueError("Discrete pure medium is incompatible with source blend.")
            local_q = np.where(is_original_b & ~is_original_a, 1.0, 0.0)
        else:
            canonical_a = min(med_a, med_b)
            canonical_b = max(med_a, med_b)
            if np.any(local_a != canonical_a) or np.any(local_b != canonical_b):
                raise ValueError("Discrete mixture pair is incompatible with source blend.")
            local_q = np.full(
                local_a.shape,
                medium.fraction_b if med_a <= med_b else 1.0 - medium.fraction_b,
                dtype=np.float64,
            )

        quantized[voxels] = local_q

    return quantized


def evaluate_discretization(
    blends: MaterialBlendArrays,
    mass_density_g_cm3: np.ndarray,
    *,
    mixture_bins: int,
    mask: np.ndarray | None = None,
) -> DiscretizationSensitivityResult:
    density = np.asarray(mass_density_g_cm3, dtype=np.float64)
    if density.shape != blends.upper_fraction.shape:
        raise ValueError("Density and material blend arrays must have identical shapes.")

    selected = np.ones(density.shape, dtype=bool) if mask is None else np.asarray(mask, dtype=bool)
    if selected.shape != density.shape:
        raise ValueError("Sensitivity mask shape does not match density field.")
    if not np.any(selected):
        raise ValueError("Sensitivity mask contains no voxels.")

    try:
        discrete = discretize_material_blends(
            blends,
            mixture_bins=mixture_bins,
        )
        pegsless = build_pegsless_media_set(discrete)
    except ValueError as exc:
        return DiscretizationSensitivityResult(
            mixture_bins=mixture_bins,
            feasible=False,
            media_count=None,
            voxel_count=int(np.count_nonzero(selected)),
            mean_abs_fraction_error=None,
            p95_abs_fraction_error=None,
            max_abs_fraction_error=None,
            rhof_min=None,
            rhof_max=None,
            rhof_p95_abs_percent=None,
            rhof_over_20_percent_fraction=None,
            error=str(exc),
        )

    _, _, original_fraction = _canonical_blend(blends)
    quantized_fraction = _quantized_fraction_field(blends, discrete)
    fraction_error = np.abs(quantized_fraction[selected] - original_fraction[selected])

    refs = np.asarray(
        [medium.bulk_density_g_cm3 for medium in pegsless.media],
        dtype=np.float64,
    )
    voxel_refs = refs[np.asarray(discrete.medium_index, dtype=np.int32) - 1]
    rhof = density[selected] / voxel_refs[selected]
    abs_percent = np.abs(rhof - 1.0) * 100.0

    return DiscretizationSensitivityResult(
        mixture_bins=mixture_bins,
        feasible=True,
        media_count=len(discrete.media),
        voxel_count=int(np.count_nonzero(selected)),
        mean_abs_fraction_error=float(np.mean(fraction_error)),
        p95_abs_fraction_error=float(np.percentile(fraction_error, 95)),
        max_abs_fraction_error=float(np.max(fraction_error)),
        rhof_min=float(np.min(rhof)),
        rhof_max=float(np.max(rhof)),
        rhof_p95_abs_percent=float(np.percentile(abs_percent, 95)),
        rhof_over_20_percent_fraction=float(np.mean(abs_percent > 20.0)),
        error=None,
    )


def run_discretization_sensitivity(
    blends: MaterialBlendArrays,
    mass_density_g_cm3: np.ndarray,
    *,
    candidate_mixture_bins: tuple[int, ...] = (5, 9, 19),
    mask: np.ndarray | None = None,
) -> tuple[DiscretizationSensitivityResult, ...]:
    if not candidate_mixture_bins:
        raise ValueError("At least one mixture-bin candidate is required.")
    if any(value < 1 for value in candidate_mixture_bins):
        raise ValueError("All mixture-bin candidates must be >= 1.")

    return tuple(
        evaluate_discretization(
            blends,
            mass_density_g_cm3,
            mixture_bins=value,
            mask=mask,
        )
        for value in candidate_mixture_bins
    )
