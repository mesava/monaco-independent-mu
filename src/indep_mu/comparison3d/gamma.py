from __future__ import annotations

from dataclasses import dataclass
import importlib

import numpy as np

from indep_mu.dicom.rtdose import RtDose

from .dose_difference import DoseDifferenceResult


@dataclass(frozen=True)
class GammaConfig:
    """Explicit gamma-analysis configuration."""

    dose_percent_threshold: float
    distance_mm_threshold: float
    lower_percent_dose_cutoff: float
    local_gamma: bool = False
    global_normalisation_gy: float | None = None
    interp_fraction: int = 10
    max_gamma: float | None = 2.0

    def __post_init__(self) -> None:
        if self.dose_percent_threshold <= 0:
            raise ValueError("dose_percent_threshold must be positive.")
        if self.distance_mm_threshold <= 0:
            raise ValueError("distance_mm_threshold must be positive.")
        if not 0 <= self.lower_percent_dose_cutoff < 100:
            raise ValueError("lower_percent_dose_cutoff must be in [0,100).")
        if self.interp_fraction < 1:
            raise ValueError("interp_fraction must be >= 1.")
        if self.max_gamma is not None and self.max_gamma <= 1:
            raise ValueError("max_gamma must be > 1 when supplied.")

        if self.local_gamma:
            if self.global_normalisation_gy is not None:
                raise ValueError(
                    "global_normalisation_gy must be omitted for local gamma."
                )
        else:
            if (
                self.global_normalisation_gy is None
                or not np.isfinite(self.global_normalisation_gy)
                or self.global_normalisation_gy <= 0
            ):
                raise ValueError(
                    "Global gamma requires an explicit positive "
                    "global_normalisation_gy."
                )


@dataclass(frozen=True)
class GammaResult:
    gamma: np.ndarray
    analysed_voxel_count: int
    passing_voxel_count: int
    pass_rate_percent: float
    config: GammaConfig


def _reference_axes_mm(reference: RtDose) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Axes in the same z,y,x order as the RTDOSE array."""

    geometry = reference.geometry
    row_spacing, column_spacing = geometry.pixel_spacing_mm

    z = np.asarray(geometry.local_frame_offsets_mm, dtype=np.float64)
    y = np.arange(geometry.rows, dtype=np.float64) * row_spacing
    x = np.arange(geometry.columns, dtype=np.float64) * column_spacing

    return z, y, x


def gamma_pass_rate(
    gamma: np.ndarray,
    *,
    mask: np.ndarray | None = None,
) -> tuple[int, int, float]:
    values = np.asarray(gamma, dtype=np.float64)
    selected = np.isfinite(values)

    if mask is not None:
        analysis_mask = np.asarray(mask, dtype=bool)
        if analysis_mask.shape != values.shape:
            raise ValueError("Gamma analysis mask shape does not match gamma grid.")
        selected &= analysis_mask

    analysed = int(np.count_nonzero(selected))
    if analysed == 0:
        raise ValueError("Gamma result contains no analysed voxels.")

    passed = int(np.count_nonzero(values[selected] <= 1.0))
    return analysed, passed, 100.0 * passed / analysed


def run_gamma_on_rtdose_grid(
    reference: RtDose,
    comparison: DoseDifferenceResult,
    config: GammaConfig,
    *,
    analysis_mask: np.ndarray | None = None,
) -> GammaResult:
    """Run PyMedPhys gamma after MC has been sampled onto the RTDOSE grid.

    The TPS RTDOSE remains the reference.  Evaluation dose is the independently
    calculated MC field sampled at the same physical reference-grid centres.

    This wrapper deliberately requires an explicit global normalisation for
    global gamma so that the project never silently falls back to reference-max
    normalisation.
    """

    if comparison.reference_dose_gy.shape != reference.shape:
        raise ValueError("Comparison/reference grid shapes do not match.")
    if comparison.evaluated_dose_gy.shape != reference.shape:
        raise ValueError("Evaluation/reference grid shapes do not match.")

    try:
        pymedphys = importlib.import_module("pymedphys")
    except ImportError as exc:
        raise RuntimeError(
            "Gamma analysis requires the optional PyMedPhys dependency."
        ) from exc

    axes = _reference_axes_mm(reference)

    gamma = np.asarray(
        pymedphys.gamma(
            axes,
            comparison.reference_dose_gy,
            axes,
            comparison.evaluated_dose_gy,
            dose_percent_threshold=config.dose_percent_threshold,
            distance_mm_threshold=config.distance_mm_threshold,
            lower_percent_dose_cutoff=config.lower_percent_dose_cutoff,
            interp_fraction=config.interp_fraction,
            max_gamma=config.max_gamma,
            local_gamma=config.local_gamma,
            global_normalisation=config.global_normalisation_gy,
            skip_once_passed=True,
        ),
        dtype=np.float64,
    )

    if gamma.shape != reference.shape:
        raise ValueError(
            f"Gamma output shape {gamma.shape} does not match RTDOSE "
            f"shape {reference.shape}."
        )

    analysed, passed, pass_rate = gamma_pass_rate(
        gamma,
        mask=analysis_mask,
    )

    return GammaResult(
        gamma=gamma,
        analysed_voxel_count=analysed,
        passing_voxel_count=passed,
        pass_rate_percent=pass_rate,
        config=config,
    )
