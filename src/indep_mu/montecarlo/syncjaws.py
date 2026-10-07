from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from indep_mu.beam_model.iec_coordinates import (
    IsocenterProjection,
    split_dicom_banks,
)
from indep_mu.dicom.rtplan import Beam


@dataclass(frozen=True)
class FocusedJawPairGeometry:
    """One source-focused SYNCJAWS pair."""

    device_type: str
    zmin_cm: float
    zmax_cm: float
    negative_bank: int

    def __post_init__(self) -> None:
        if self.zmin_cm <= 0 or self.zmax_cm <= self.zmin_cm:
            raise ValueError("Jaw geometry requires 0 < zmin_cm < zmax_cm.")
        if self.negative_bank not in {1, 2}:
            raise ValueError("negative_bank must be 1 or 2.")


@dataclass(frozen=True)
class SyncJawsPairState:
    zmin_cm: float
    zmax_cm: float
    positive_front_cm: float
    positive_back_cm: float
    negative_front_cm: float
    negative_back_cm: float


@dataclass(frozen=True)
class SyncJawsSequencePoint:
    mu_index: float
    pairs: tuple[SyncJawsPairState, ...]


@dataclass(frozen=True)
class SyncJawsSequence:
    points: tuple[SyncJawsSequencePoint, ...]
    geometries: tuple[FocusedJawPairGeometry, ...]


def _pair_state(
    *,
    positions_mm: tuple[float, ...],
    geometry: FocusedJawPairGeometry,
    sad_mm: float,
) -> SyncJawsPairState:
    banks = split_dicom_banks(positions_mm, pair_count=1)
    bank1 = float(banks.bank1_mm[0])
    bank2 = float(banks.bank2_mm[0])

    if geometry.negative_bank == 1:
        negative_iso, positive_iso = bank1, bank2
    else:
        negative_iso, positive_iso = bank2, bank1

    projection = IsocenterProjection(sad_mm)

    return SyncJawsPairState(
        zmin_cm=geometry.zmin_cm,
        zmax_cm=geometry.zmax_cm,
        positive_front_cm=float(
            projection.to_source_distance_cm(
                positive_iso,
                source_distance_cm=geometry.zmin_cm,
            )
        ),
        positive_back_cm=float(
            projection.to_source_distance_cm(
                positive_iso,
                source_distance_cm=geometry.zmax_cm,
            )
        ),
        negative_front_cm=float(
            projection.to_source_distance_cm(
                negative_iso,
                source_distance_cm=geometry.zmin_cm,
            )
        ),
        negative_back_cm=float(
            projection.to_source_distance_cm(
                negative_iso,
                source_distance_cm=geometry.zmax_cm,
            )
        ),
    )


def build_syncjaws_sequence(
    beam: Beam,
    *,
    geometries: tuple[FocusedJawPairGeometry, ...],
    equal_mu_tolerance: float = 1e-10,
) -> SyncJawsSequence:
    """Build a source-focused SYNCJAWS sequence for one beam."""

    if beam.source_axis_distance_mm is None:
        raise ValueError("Beam SourceAxisDistance is required.")
    if not geometries:
        raise ValueError("At least one jaw geometry is required.")

    points: list[SyncJawsSequencePoint] = []
    for cp in beam.control_points:
        pair_states = tuple(
            _pair_state(
                positions_mm=cp.positions(geometry.device_type),
                geometry=geometry,
                sad_mm=beam.source_axis_distance_mm,
            )
            for geometry in geometries
        )
        point = SyncJawsSequencePoint(
            mu_index=(
                cp.cumulative_meterset_weight
                / beam.final_cumulative_meterset_weight
            ),
            pairs=pair_states,
        )
        if points and abs(point.mu_index - points[-1].mu_index) <= equal_mu_tolerance:
            points[-1] = point
        else:
            points.append(point)

    indices = np.asarray([point.mu_index for point in points], dtype=np.float64)
    if len(points) < 2:
        raise ValueError("SYNCJAWS dynamic sequence requires at least two states.")
    if abs(indices[0]) > equal_mu_tolerance:
        raise ValueError("SYNCJAWS sequence must start at index 0.")
    if abs(indices[-1] - 1.0) > equal_mu_tolerance:
        raise ValueError("SYNCJAWS sequence must end at index 1.")
    if np.any(np.diff(indices) <= equal_mu_tolerance):
        raise ValueError("SYNCJAWS indices must be strictly increasing.")

    return SyncJawsSequence(
        points=tuple(points),
        geometries=geometries,
    )


def render_syncjaws_sequence(sequence: SyncJawsSequence) -> str:
    """Render the external opening file consumed by BEAMnrc SYNCJAWS."""

    lines = [f"{len(sequence.points):10d}"]
    for point in sequence.points:
        lines.append(f"{point.mu_index:15.8f}")
        for pair in point.pairs:
            lines.append(
                f"{pair.zmin_cm:15.8f}"
                f"{pair.zmax_cm:15.8f}"
                f"{pair.positive_front_cm:15.8f}"
                f"{pair.positive_back_cm:15.8f}"
                f"{pair.negative_front_cm:15.8f}"
                f"{pair.negative_back_cm:15.8f}"
            )
    return "\n".join(lines) + "\n"


def write_syncjaws_sequence(path: str | Path, sequence: SyncJawsSequence) -> None:
    Path(path).write_text(
        render_syncjaws_sequence(sequence),
        encoding="ascii",
    )
