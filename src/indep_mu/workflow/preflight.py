from __future__ import annotations

from dataclasses import dataclass

from indep_mu.dicom.ct import CtSeries
from indep_mu.dicom.rtplan import RtPlan
from indep_mu.montecarlo.source21_orientation import (
    ZhanHfsSource21OrientationMapper,
)


@dataclass(frozen=True)
class TransportPreflightMessage:
    severity: str
    code: str
    message: str
    beam_number: int | None = None
    control_point_index: int | None = None


@dataclass(frozen=True)
class TransportPreflight:
    messages: tuple[TransportPreflightMessage, ...]

    @property
    def errors(self) -> tuple[TransportPreflightMessage, ...]:
        return tuple(item for item in self.messages if item.severity == "ERROR")

    @property
    def warnings(self) -> tuple[TransportPreflightMessage, ...]:
        return tuple(item for item in self.messages if item.severity == "WARNING")

    @property
    def passed(self) -> bool:
        return not self.errors


def run_transport_preflight(
    ct: CtSeries,
    plan: RtPlan,
    *,
    calibration_kvp: float = 120.0,
    kvp_tolerance: float = 0.5,
) -> TransportPreflight:
    """Check case features that gate the current HFS/source21 transport path."""

    messages: list[TransportPreflightMessage] = []

    patient_position = ct.geometry.patient_position
    if patient_position is None:
        messages.append(
            TransportPreflightMessage(
                "ERROR",
                "PATIENT_POSITION_MISSING",
                "CT PatientPosition is missing; current source21 mapping cannot "
                "prove the patient-to-machine orientation.",
            )
        )
    elif patient_position.upper() != "HFS":
        messages.append(
            TransportPreflightMessage(
                "ERROR",
                "PATIENT_POSITION_UNSUPPORTED",
                f"PatientPosition={patient_position!r}; current source21 mapper "
                "is validated only for HFS.",
            )
        )

    if ct.geometry.kvp is None:
        messages.append(
            TransportPreflightMessage(
                "ERROR",
                "CT_KVP_MISSING",
                "CT KVP is missing; the 120 kV scanner calibration cannot be "
                "selected safely.",
            )
        )
    elif abs(ct.geometry.kvp - calibration_kvp) > kvp_tolerance:
        messages.append(
            TransportPreflightMessage(
                "ERROR",
                "CT_CALIBRATION_KVP_MISMATCH",
                f"CT KVP={ct.geometry.kvp:g} kV does not match the configured "
                f"{calibration_kvp:g} kV calibration.",
            )
        )

    if patient_position is not None and patient_position.upper() == "HFS":
        mapper = ZhanHfsSource21OrientationMapper(patient_position="HFS")
        for beam in plan.beams:
            for cp in beam.control_points:
                if (
                    cp.gantry_angle_deg is None
                    or cp.collimator_angle_deg is None
                    or cp.patient_support_angle_deg is None
                ):
                    # Missing angles are already handled by RTPLAN preflight.
                    continue
                try:
                    mapper.map_angles(
                        gantry_deg=cp.gantry_angle_deg,
                        collimator_deg=cp.collimator_angle_deg,
                        patient_support_deg=cp.patient_support_angle_deg,
                    )
                except ValueError as exc:
                    messages.append(
                        TransportPreflightMessage(
                            "ERROR",
                            "SOURCE21_ORIENTATION_UNSUPPORTED",
                            str(exc),
                            beam.number,
                            cp.index,
                        )
                    )

    return TransportPreflight(messages=tuple(messages))
