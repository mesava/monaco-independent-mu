from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from indep_mu.dicom.rtplan import Beam


@dataclass(frozen=True)
class GeometryMessage:
    severity: str
    code: str
    message: str


@dataclass(frozen=True)
class AgilityDicomGeometryReport:
    messages: tuple[GeometryMessage, ...]
    mlc_device_type: str | None
    leaf_pairs: int | None
    field_span_mm: float | None
    nominal_leaf_width_mm: float | None

    @property
    def errors(self) -> tuple[GeometryMessage, ...]:
        return tuple(item for item in self.messages if item.severity == "ERROR")

    @property
    def passed(self) -> bool:
        return not self.errors


def validate_agility_dicom_geometry(
    beam: Beam,
    *,
    expected_leaf_pairs: int = 80,
    expected_leaf_width_mm: float = 5.0,
    expected_field_span_mm: float = 400.0,
    tolerance_mm: float = 1e-3,
) -> AgilityDicomGeometryReport:
    """Validate delivery-coordinate geometry exported for an Agility MLC.

    This validates DICOM delivery geometry only. It does NOT claim that DICOM
    leaf boundaries or Monaco tuning parameters define the physical MC leaf
    body, rounded tip, tongue-and-groove or material composition.
    """

    messages: list[GeometryMessage] = []

    mlcs = [
        item for item in beam.device_definitions
        if item.device_type.startswith("MLC")
    ]
    if len(mlcs) != 1:
        messages.append(
            GeometryMessage(
                "ERROR",
                "MLC_DEVICE_COUNT",
                f"Expected exactly one MLC device, found {len(mlcs)}.",
            )
        )
        return AgilityDicomGeometryReport(
            messages=tuple(messages),
            mlc_device_type=None,
            leaf_pairs=None,
            field_span_mm=None,
            nominal_leaf_width_mm=None,
        )

    mlc = mlcs[0]
    if mlc.number_of_leaf_jaw_pairs != expected_leaf_pairs:
        messages.append(
            GeometryMessage(
                "ERROR",
                "LEAF_PAIR_COUNT",
                f"Expected {expected_leaf_pairs} Agility leaf pairs, "
                f"DICOM contains {mlc.number_of_leaf_jaw_pairs}.",
            )
        )

    boundaries = mlc.leaf_position_boundaries_mm
    if boundaries is None:
        messages.append(
            GeometryMessage(
                "ERROR",
                "MISSING_LEAF_BOUNDARIES",
                "MLC LeafPositionBoundaries are missing.",
            )
        )
        return AgilityDicomGeometryReport(
            messages=tuple(messages),
            mlc_device_type=mlc.device_type,
            leaf_pairs=mlc.number_of_leaf_jaw_pairs,
            field_span_mm=None,
            nominal_leaf_width_mm=None,
        )

    boundary_values = np.asarray(boundaries, dtype=np.float64)
    widths = np.diff(boundary_values)
    span = float(boundary_values[-1] - boundary_values[0])
    nominal_width = float(np.median(widths))

    if np.any(widths <= 0):
        messages.append(
            GeometryMessage(
                "ERROR",
                "NON_MONOTONIC_LEAF_BOUNDARIES",
                "LeafPositionBoundaries are not strictly increasing.",
            )
        )
    if not np.allclose(
        widths,
        expected_leaf_width_mm,
        rtol=0.0,
        atol=tolerance_mm,
    ):
        messages.append(
            GeometryMessage(
                "ERROR",
                "LEAF_WIDTH",
                f"Agility expected {expected_leaf_width_mm} mm projected width; "
                f"observed range {float(widths.min()):.6g}–"
                f"{float(widths.max()):.6g} mm.",
            )
        )
    if abs(span - expected_field_span_mm) > tolerance_mm:
        messages.append(
            GeometryMessage(
                "ERROR",
                "MLC_FIELD_SPAN",
                f"Expected field span {expected_field_span_mm} mm, got {span} mm.",
            )
        )

    if beam.source_axis_distance_mm is not None and abs(
        beam.source_axis_distance_mm - 1000.0
    ) > tolerance_mm:
        messages.append(
            GeometryMessage(
                "ERROR",
                "SAD",
                f"Expected Versa HD SAD 1000 mm, got "
                f"{beam.source_axis_distance_mm} mm.",
            )
        )

    if not any(
        item.device_type in {"ASYMY", "Y"}
        for item in beam.device_definitions
    ):
        messages.append(
            GeometryMessage(
                "WARNING",
                "Y_JAWS_NOT_IDENTIFIED",
                "No Y-jaw device (ASYMY/Y) identified in DICOM.",
            )
        )

    return AgilityDicomGeometryReport(
        messages=tuple(messages),
        mlc_device_type=mlc.device_type,
        leaf_pairs=mlc.number_of_leaf_jaw_pairs,
        field_span_mm=span,
        nominal_leaf_width_mm=nominal_width,
    )
