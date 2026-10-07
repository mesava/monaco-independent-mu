from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from indep_mu.dicom.rtplan import Beam, BeamLimitingDeviceState


@dataclass(frozen=True)
class DeliveryInterpolationPolicy:
    """Explicit non-DICOM assumption for states between control points.

    DICOM defines the states at control points and rotation direction for the
    following segment, but does not define machine behaviour between control
    points. Therefore there is intentionally no implicit/default policy.
    """

    policy_id: str
    evidence: str
    linear_in_cumulative_meterset: bool = True

    def __post_init__(self) -> None:
        if not self.policy_id.strip():
            raise ValueError("policy_id is required.")
        if not self.evidence.strip():
            raise ValueError(
                "Interpolation policy requires an evidence/validation statement."
            )
        if not self.linear_in_cumulative_meterset:
            raise ValueError(
                "Only linear-in-CMW policy is implemented; other policies "
                "require a separate implementation."
            )


@dataclass(frozen=True)
class SampledBeamState:
    segment_index: int
    fraction: float
    cumulative_meterset_weight: float
    mu_from_beam_start: float
    gantry_angle_deg: float | None
    collimator_angle_deg: float | None
    patient_support_angle_deg: float | None
    device_positions: tuple[BeamLimitingDeviceState, ...]

    def positions(self, device_type: str) -> tuple[float, ...]:
        for item in self.device_positions:
            if item.device_type == device_type:
                return item.positions_mm
        raise KeyError(device_type)


def directed_rotation_delta_deg(
    start_deg: float,
    end_deg: float,
    direction: str | None,
    *,
    tolerance_deg: float = 1e-9,
) -> float:
    """Return signed DICOM rotation from one CP to the next.

    For machine rotation attributes, DICOM defines CC as the direction of
    increasing IEC angle and CW as decreasing. Equal start/end angles with CW
    or CC encode a full 360-degree rotation; NONE encodes no movement.
    """

    start = float(start_deg) % 360.0
    end = float(end_deg) % 360.0
    mode = "NONE" if direction is None else direction.upper()

    if mode == "NONE":
        raw = abs(end - start)
        circular = min(raw, 360.0 - raw)
        if circular > tolerance_deg:
            raise ValueError(
                f"RotationDirection=NONE but angle changes {start} -> {end} deg."
            )
        return 0.0

    if mode == "CC":
        delta = (end - start) % 360.0
        if delta <= tolerance_deg:
            delta = 360.0
        return delta

    if mode == "CW":
        magnitude = (start - end) % 360.0
        if magnitude <= tolerance_deg:
            magnitude = 360.0
        return -magnitude

    raise ValueError(f"Unsupported DICOM rotation direction {direction!r}.")


def _interpolate_angle(
    start: float | None,
    end: float | None,
    direction: str | None,
    fraction: float,
    name: str,
) -> float | None:
    if start is None and end is None:
        return None
    if start is None or end is None:
        raise ValueError(f"{name} is not available at both segment endpoints.")
    delta = directed_rotation_delta_deg(start, end, direction)
    return float((start + fraction * delta) % 360.0)


def sample_beam_segment(
    beam: Beam,
    segment_index: int,
    fraction: float,
    *,
    policy: DeliveryInterpolationPolicy,
) -> SampledBeamState:
    """Sample one dynamic interval using an explicitly selected policy.

    Linear interpolation in cumulative meterset is useful for sequence-file
    generation and numerical quadrature, but it is a machine-model assumption,
    not a rule guaranteed by DICOM itself.
    """

    if policy is None:
        raise ValueError("An explicit DeliveryInterpolationPolicy is required.")
    if not 0.0 <= fraction <= 1.0:
        raise ValueError("fraction must be in [0, 1].")

    segments = beam.segments
    if not 0 <= segment_index < len(segments):
        raise IndexError("segment_index out of range.")

    start = beam.control_points[segment_index]
    end = beam.control_points[segment_index + 1]
    segment = segments[segment_index]

    states: list[BeamLimitingDeviceState] = []
    for definition in beam.device_definitions:
        device_type = definition.device_type
        start_values = np.asarray(start.positions(device_type), dtype=np.float64)
        end_values = np.asarray(end.positions(device_type), dtype=np.float64)
        if start_values.shape != end_values.shape:
            raise ValueError(
                f"{device_type} endpoint position arrays have different shapes."
            )
        values = start_values + fraction * (end_values - start_values)
        states.append(
            BeamLimitingDeviceState(
                device_type=device_type,
                positions_mm=tuple(float(value) for value in values),
            )
        )

    cmw = (
        start.cumulative_meterset_weight
        + fraction * segment.delta_cumulative_meterset_weight
    )
    mu = beam.beam_meterset_mu * cmw / beam.final_cumulative_meterset_weight

    return SampledBeamState(
        segment_index=segment_index,
        fraction=float(fraction),
        cumulative_meterset_weight=float(cmw),
        mu_from_beam_start=float(mu),
        gantry_angle_deg=_interpolate_angle(
            start.gantry_angle_deg,
            end.gantry_angle_deg,
            start.gantry_rotation_direction,
            fraction,
            "GantryAngle",
        ),
        collimator_angle_deg=_interpolate_angle(
            start.collimator_angle_deg,
            end.collimator_angle_deg,
            start.collimator_rotation_direction,
            fraction,
            "BeamLimitingDeviceAngle",
        ),
        patient_support_angle_deg=_interpolate_angle(
            start.patient_support_angle_deg,
            end.patient_support_angle_deg,
            start.patient_support_rotation_direction,
            fraction,
            "PatientSupportAngle",
        ),
        device_positions=tuple(states),
    )
