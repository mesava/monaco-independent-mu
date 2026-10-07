from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Protocol, runtime_checkable

import numpy as np

from indep_mu.beam_model.agility_rounded_tip import RoundedLeafTipTangentGeometry
from indep_mu.beam_model.iec_coordinates import (
    DicomBankPositions,
    IsocenterProjection,
    split_dicom_banks,
)
from indep_mu.dicom.rtplan import Beam


@dataclass(frozen=True)
class SyncMlceOpening:
    negative_cm: np.ndarray
    positive_cm: np.ndarray

    def __post_init__(self) -> None:
        negative = np.asarray(self.negative_cm, dtype=np.float64)
        positive = np.asarray(self.positive_cm, dtype=np.float64)
        if negative.ndim != 1 or positive.ndim != 1:
            raise ValueError("SYNCMLCE opening arrays must be one-dimensional.")
        if negative.shape != positive.shape:
            raise ValueError("SYNCMLCE opening arrays must have equal length.")
        if not np.all(np.isfinite(negative)) or not np.all(np.isfinite(positive)):
            raise ValueError("SYNCMLCE opening coordinates must be finite.")


@runtime_checkable
class SyncMlceCoordinateMapper(Protocol):
    """Map DICOM bank positions to SYNCMLCE NEG/POS coordinates."""

    def map_banks(self, banks: DicomBankPositions) -> SyncMlceOpening:
        ...


@dataclass(frozen=True)
class SourceFocusedSyncMlceMapper:
    """Exact mapper for an ENDTYPE=1 leaf end focused to the source.

    SYNCMLCE ENDTYPE=1 expects the opening coordinate at ZMIN.  DICOM
    isocenter-projected positions can therefore be scaled by ZMIN/SAD.

    This mapper must NOT be used for rounded/cylindrical Agility leaf tips.
    ENDTYPE=0 expects cylinder-origin coordinates instead of a projected front
    opening and requires a separately validated physical leaf-tip model.
    """

    sad_mm: float
    zmin_cm: float
    negative_bank: int

    def __post_init__(self) -> None:
        if self.zmin_cm <= 0:
            raise ValueError("zmin_cm must be positive.")
        if self.negative_bank not in {1, 2}:
            raise ValueError("negative_bank must be 1 or 2.")

    def map_banks(self, banks: DicomBankPositions) -> SyncMlceOpening:
        projection = IsocenterProjection(self.sad_mm)
        bank1 = projection.to_source_distance_cm(
            banks.bank1_mm, source_distance_cm=self.zmin_cm
        )
        bank2 = projection.to_source_distance_cm(
            banks.bank2_mm, source_distance_cm=self.zmin_cm
        )
        if self.negative_bank == 1:
            negative, positive = bank1, bank2
        else:
            negative, positive = bank2, bank1
        return SyncMlceOpening(
            negative_cm=np.asarray(negative, dtype=np.float64),
            positive_cm=np.asarray(positive, dtype=np.float64),
        )



@dataclass(frozen=True)
class ResearchRoundedTipSyncMlceMapper:
    """Research-only DICOM edge -> SYNCMLCE ENDTYPE=0 cylinder mapper.

    This uses exact source-ray tangency for a cylindrical tip.  It is not a
    commissioned Agility mapper and must not be selected by a clinical default.

    Bank identity is explicit because DICOM bank order must never be inferred
    from the current sign of a moving/interdigitating leaf position.
    """

    sad_mm: float
    radius_cm: float
    cylinder_axis_z_cm: float
    negative_bank: int

    def __post_init__(self) -> None:
        if not np.isfinite(self.sad_mm) or self.sad_mm <= 0:
            raise ValueError("sad_mm must be finite and positive.")
        if self.negative_bank not in {1, 2}:
            raise ValueError("negative_bank must be 1 or 2.")

        # Validate physical radius/CIL constraints immediately.
        RoundedLeafTipTangentGeometry(
            sad_cm=self.sad_mm / 10.0,
            radius_cm=self.radius_cm,
            cylinder_axis_z_cm=self.cylinder_axis_z_cm,
        )

    def map_banks(self, banks: DicomBankPositions) -> SyncMlceOpening:
        geometry = RoundedLeafTipTangentGeometry(
            sad_cm=self.sad_mm / 10.0,
            radius_cm=self.radius_cm,
            cylinder_axis_z_cm=self.cylinder_axis_z_cm,
        )

        bank1_iso_cm = np.asarray(banks.bank1_mm, dtype=np.float64) / 10.0
        bank2_iso_cm = np.asarray(banks.bank2_mm, dtype=np.float64) / 10.0

        if self.negative_bank == 1:
            negative_edge = bank1_iso_cm
            positive_edge = bank2_iso_cm
        else:
            negative_edge = bank2_iso_cm
            positive_edge = bank1_iso_cm

        negative_origin = geometry.projected_edge_to_cylinder_origin_cm(
            negative_edge,
            opening_side="negative",
        )
        positive_origin = geometry.projected_edge_to_cylinder_origin_cm(
            positive_edge,
            opening_side="positive",
        )

        return SyncMlceOpening(
            negative_cm=np.asarray(negative_origin, dtype=np.float64),
            positive_cm=np.asarray(positive_origin, dtype=np.float64),
        )


@dataclass(frozen=True)
class SyncMlceSequencePoint:
    mu_index: float
    opening: SyncMlceOpening


@dataclass(frozen=True)
class SyncMlceSequence:
    title: str
    device_type: str
    leaf_pairs: int
    points: tuple[SyncMlceSequencePoint, ...]


def _coalesce_equal_mu_points(
    points: list[SyncMlceSequencePoint],
    *,
    tolerance: float,
) -> tuple[SyncMlceSequencePoint, ...]:
    """Keep the last state at a zero-MU transition.

    Equal CMW control points represent no delivered dose between the states.
    SYNCMLCE dynamic interpolation requires a usable ordered index; retaining
    the last state makes the instantaneous zero-weight transition explicit.
    """

    result: list[SyncMlceSequencePoint] = []
    for point in points:
        if result and abs(point.mu_index - result[-1].mu_index) <= tolerance:
            result[-1] = point
        else:
            result.append(point)
    return tuple(result)


def build_syncmlce_sequence(
    beam: Beam,
    *,
    mapper: SyncMlceCoordinateMapper,
    device_type: str | None = None,
    equal_mu_tolerance: float = 1e-10,
) -> SyncMlceSequence:
    """Build a per-beam SYNCMLCE sequence from DICOM control points.

    MUINDEX is normalized within this beam: CMW / FinalCMW.  This matches the
    project's per-beam Monte Carlo strategy; absolute beam MU is handled by
    dose normalization and plan summation, not by changing SYNCMLCE MUINDEX.
    """

    mlc_definitions = [
        item
        for item in beam.device_definitions
        if item.device_type.startswith("MLC")
    ]
    if device_type is None:
        if len(mlc_definitions) != 1:
            raise ValueError(
                f"Expected exactly one MLC device, found {len(mlc_definitions)}."
            )
        definition = mlc_definitions[0]
    else:
        definition = beam.device_definition(device_type)
        if not definition.device_type.startswith("MLC"):
            raise ValueError(f"{device_type} is not an MLC device.")

    if beam.final_cumulative_meterset_weight <= 0:
        raise ValueError("FinalCumulativeMetersetWeight must be positive.")

    points: list[SyncMlceSequencePoint] = []
    for cp in beam.control_points:
        banks = split_dicom_banks(
            cp.positions(definition.device_type),
            pair_count=definition.number_of_leaf_jaw_pairs,
        )
        opening = mapper.map_banks(banks)
        if opening.negative_cm.size != definition.number_of_leaf_jaw_pairs:
            raise ValueError("Coordinate mapper returned wrong number of leaf pairs.")

        points.append(
            SyncMlceSequencePoint(
                mu_index=(
                    cp.cumulative_meterset_weight
                    / beam.final_cumulative_meterset_weight
                ),
                opening=opening,
            )
        )

    sequence_points = _coalesce_equal_mu_points(
        points,
        tolerance=equal_mu_tolerance,
    )
    indices = np.asarray([point.mu_index for point in sequence_points])
    if indices.size < 2:
        raise ValueError("SYNCMLCE dynamic sequence requires at least two MU states.")
    if abs(indices[0]) > equal_mu_tolerance:
        raise ValueError("SYNCMLCE sequence must start at MUINDEX 0.")
    if abs(indices[-1] - 1.0) > equal_mu_tolerance:
        raise ValueError("SYNCMLCE sequence must end at MUINDEX 1.")
    if np.any(np.diff(indices) <= equal_mu_tolerance):
        raise ValueError("SYNCMLCE MUINDEX values must be strictly increasing.")

    return SyncMlceSequence(
        title=f"Beam {beam.number} {beam.name}"[:80],
        device_type=definition.device_type,
        leaf_pairs=definition.number_of_leaf_jaw_pairs,
        points=sequence_points,
    )


def render_syncmlce_sequence(sequence: SyncMlceSequence) -> str:
    """Render the external leaf-opening file consumed by BEAMnrc SYNCMLCE."""

    lines = [sequence.title[:80], f"{len(sequence.points):10d}"]
    for point in sequence.points:
        lines.append(f"{point.mu_index:15.8f}")
        for negative, positive in zip(
            point.opening.negative_cm,
            point.opening.positive_cm,
        ):
            lines.append(f"{negative:15.8f}{positive:15.8f}{1:5d}")
    return "\n".join(lines) + "\n"


def write_syncmlce_sequence(
    path: str | Path,
    sequence: SyncMlceSequence,
) -> None:
    Path(path).write_text(
        render_syncmlce_sequence(sequence),
        encoding="ascii",
    )
