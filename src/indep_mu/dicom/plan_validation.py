from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .rtplan import Beam, RtPlan


@dataclass(frozen=True)
class BeamDeliveryFeatures:
    beam_number: int
    delivery_class: str
    control_point_count: int
    segment_count: int
    beam_meterset_mu: float
    gantry_dynamic: bool
    collimator_dynamic: bool
    couch_dynamic: bool
    jaw_dynamic: bool
    mlc_dynamic: bool
    energies_mv: tuple[float, ...]
    device_types: tuple[str, ...]
    isocenter_count: int


@dataclass(frozen=True)
class PreflightMessage:
    severity: str
    code: str
    message: str
    beam_number: int | None = None


@dataclass(frozen=True)
class PlanPreflight:
    beam_features: tuple[BeamDeliveryFeatures, ...]
    messages: tuple[PreflightMessage, ...]

    @property
    def errors(self) -> tuple[PreflightMessage, ...]:
        return tuple(item for item in self.messages if item.severity == "ERROR")

    @property
    def warnings(self) -> tuple[PreflightMessage, ...]:
        return tuple(item for item in self.messages if item.severity == "WARNING")

    @property
    def passed(self) -> bool:
        return not self.errors


def _circular_change(values: list[float], tolerance_deg: float) -> bool:
    if len(values) < 2:
        return False
    first = float(values[0]) % 360.0
    for value in values[1:]:
        current = float(value) % 360.0
        delta = abs(current - first)
        delta = min(delta, 360.0 - delta)
        if delta > tolerance_deg:
            return True
    return False


def _linear_change(values: list[float], tolerance: float) -> bool:
    if len(values) < 2:
        return False
    return bool(np.ptp(np.asarray(values, dtype=np.float64)) > tolerance)


def _device_dynamic(
    beam: Beam,
    *,
    mlc: bool,
    tolerance_mm: float,
) -> bool:
    target_types = [
        item.device_type
        for item in beam.device_definitions
        if item.device_type.startswith("MLC") == mlc
    ]
    for device_type in target_types:
        reference = np.asarray(
            beam.control_points[0].positions(device_type),
            dtype=np.float64,
        )
        for cp in beam.control_points[1:]:
            current = np.asarray(cp.positions(device_type), dtype=np.float64)
            if not np.allclose(
                current,
                reference,
                rtol=0.0,
                atol=tolerance_mm,
            ):
                return True
    return False


def analyze_beam_delivery(
    beam: Beam,
    *,
    angle_tolerance_deg: float = 1e-4,
    position_tolerance_mm: float = 1e-4,
) -> BeamDeliveryFeatures:
    gantry_values = [
        cp.gantry_angle_deg
        for cp in beam.control_points
        if cp.gantry_angle_deg is not None
    ]
    collimator_values = [
        cp.collimator_angle_deg
        for cp in beam.control_points
        if cp.collimator_angle_deg is not None
    ]
    couch_values = [
        cp.patient_support_angle_deg
        for cp in beam.control_points
        if cp.patient_support_angle_deg is not None
    ]

    gantry_dynamic = _circular_change(gantry_values, angle_tolerance_deg)
    collimator_dynamic = _circular_change(
        collimator_values,
        angle_tolerance_deg,
    )
    couch_dynamic = _circular_change(couch_values, angle_tolerance_deg)
    mlc_dynamic = _device_dynamic(
        beam,
        mlc=True,
        tolerance_mm=position_tolerance_mm,
    )
    jaw_dynamic = _device_dynamic(
        beam,
        mlc=False,
        tolerance_mm=position_tolerance_mm,
    )

    if gantry_dynamic and mlc_dynamic:
        delivery_class = "VMAT"
    elif gantry_dynamic:
        delivery_class = "ARC"
    elif mlc_dynamic:
        delivery_class = "DYNAMIC_MLC"
    elif jaw_dynamic:
        delivery_class = "DYNAMIC_JAW"
    else:
        delivery_class = "STATIC"

    energies = sorted(
        {
            float(cp.nominal_beam_energy_mv)
            for cp in beam.control_points
            if cp.nominal_beam_energy_mv is not None
        }
    )
    isocenters = {
        tuple(round(value, 6) for value in cp.isocenter_position_mm)
        for cp in beam.control_points
        if cp.isocenter_position_mm is not None
    }

    return BeamDeliveryFeatures(
        beam_number=beam.number,
        delivery_class=delivery_class,
        control_point_count=len(beam.control_points),
        segment_count=len(beam.segments),
        beam_meterset_mu=beam.beam_meterset_mu,
        gantry_dynamic=gantry_dynamic,
        collimator_dynamic=collimator_dynamic,
        couch_dynamic=couch_dynamic,
        jaw_dynamic=jaw_dynamic,
        mlc_dynamic=mlc_dynamic,
        energies_mv=tuple(energies),
        device_types=tuple(
            item.device_type for item in beam.device_definitions
        ),
        isocenter_count=len(isocenters),
    )


def run_plan_preflight(plan: RtPlan) -> PlanPreflight:
    features: list[BeamDeliveryFeatures] = []
    messages: list[PreflightMessage] = []

    for beam in plan.beams:
        item = analyze_beam_delivery(beam)
        features.append(item)

        if beam.radiation_type != "PHOTON":
            messages.append(
                PreflightMessage(
                    "ERROR",
                    "UNSUPPORTED_RADIATION_TYPE",
                    f"RadiationType={beam.radiation_type!r} is not supported.",
                    beam.number,
                )
            )

        if beam.treatment_delivery_type != "TREATMENT":
            messages.append(
                PreflightMessage(
                    "WARNING",
                    "NON_TREATMENT_BEAM",
                    "Beam is not marked as TREATMENT.",
                    beam.number,
                )
            )

        modifiers = {
            "wedges": beam.number_of_wedges,
            "compensators": beam.number_of_compensators,
            "boli": beam.number_of_boli,
            "blocks": beam.number_of_blocks,
        }
        active_modifiers = [
            f"{name}={count}"
            for name, count in modifiers.items()
            if count > 0
        ]
        if active_modifiers:
            messages.append(
                PreflightMessage(
                    "ERROR",
                    "UNSUPPORTED_BEAM_MODIFIER",
                    "Unsupported beam modifier(s): "
                    + ", ".join(active_modifiers),
                    beam.number,
                )
            )

        if len(item.energies_mv) == 0:
            messages.append(
                PreflightMessage(
                    "ERROR",
                    "MISSING_BEAM_ENERGY",
                    "No NominalBeamEnergy is available after CP inheritance.",
                    beam.number,
                )
            )
        elif len(item.energies_mv) > 1:
            messages.append(
                PreflightMessage(
                    "ERROR",
                    "ENERGY_SWITCHING",
                    f"Multiple NominalBeamEnergy values: {item.energies_mv}.",
                    beam.number,
                )
            )

        if item.isocenter_count == 0:
            messages.append(
                PreflightMessage(
                    "ERROR",
                    "MISSING_ISOCENTER",
                    "No IsocenterPosition is available.",
                    beam.number,
                )
            )
        elif item.isocenter_count > 1:
            messages.append(
                PreflightMessage(
                    "ERROR",
                    "DYNAMIC_ISOCENTER",
                    "IsocenterPosition changes within the beam.",
                    beam.number,
                )
            )

        if item.couch_dynamic:
            messages.append(
                PreflightMessage(
                    "ERROR",
                    "DYNAMIC_COUCH",
                    "Dynamic PatientSupportAngle is not supported by the current MC backend.",
                    beam.number,
                )
            )

        if item.collimator_dynamic:
            messages.append(
                PreflightMessage(
                    "ERROR",
                    "DYNAMIC_COLLIMATOR",
                    "Dynamic collimator rotation is not supported by the current MC backend.",
                    beam.number,
                )
            )

        if beam.beam_meterset_mu <= 0:
            messages.append(
                PreflightMessage(
                    "ERROR",
                    "NON_POSITIVE_MU",
                    f"BeamMeterset={beam.beam_meterset_mu}.",
                    beam.number,
                )
            )

        if not any(
            device.startswith("MLC")
            for device in item.device_types
        ):
            messages.append(
                PreflightMessage(
                    "WARNING",
                    "NO_MLC_DEVICE",
                    "Beam contains no MLC device.",
                    beam.number,
                )
            )

    return PlanPreflight(
        beam_features=tuple(features),
        messages=tuple(messages),
    )
