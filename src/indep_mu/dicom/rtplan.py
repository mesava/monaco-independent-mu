from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pydicom


@dataclass(frozen=True)
class BeamLimitingDeviceDefinition:
    device_type: str
    number_of_leaf_jaw_pairs: int
    leaf_position_boundaries_mm: tuple[float, ...] | None
    source_to_device_distance_mm: float | None


@dataclass(frozen=True)
class BeamLimitingDeviceState:
    device_type: str
    positions_mm: tuple[float, ...]


@dataclass(frozen=True)
class ControlPoint:
    index: int
    cumulative_meterset_weight: float
    gantry_angle_deg: float | None
    gantry_rotation_direction: str | None
    collimator_angle_deg: float | None
    collimator_rotation_direction: str | None
    patient_support_angle_deg: float | None
    patient_support_rotation_direction: str | None
    table_top_eccentric_angle_deg: float | None
    table_top_eccentric_rotation_direction: str | None
    table_top_vertical_position_mm: float | None
    table_top_longitudinal_position_mm: float | None
    table_top_lateral_position_mm: float | None
    nominal_beam_energy_mv: float | None
    dose_rate_set_mu_min: float | None
    source_to_surface_distance_mm: float | None
    isocenter_position_mm: tuple[float, float, float] | None
    device_positions: tuple[BeamLimitingDeviceState, ...]

    def positions(self, device_type: str) -> tuple[float, ...]:
        for item in self.device_positions:
            if item.device_type == device_type:
                return item.positions_mm
        raise KeyError(device_type)


@dataclass(frozen=True)
class DeliverySegment:
    start_control_point_index: int
    end_control_point_index: int
    delta_cumulative_meterset_weight: float
    delta_mu: float


@dataclass(frozen=True)
class Beam:
    number: int
    name: str
    description: str | None
    beam_type: str
    radiation_type: str
    treatment_delivery_type: str
    treatment_machine_name: str | None
    source_axis_distance_mm: float | None
    final_cumulative_meterset_weight: float
    beam_meterset_mu: float
    fluence_mode: str | None
    fluence_mode_id: str | None
    device_definitions: tuple[BeamLimitingDeviceDefinition, ...]
    control_points: tuple[ControlPoint, ...]

    @property
    def segments(self) -> tuple[DeliverySegment, ...]:
        result: list[DeliverySegment] = []
        final_weight = self.final_cumulative_meterset_weight
        if final_weight <= 0:
            raise ValueError(
                f"Beam {self.number} has non-positive FinalCumulativeMetersetWeight."
            )

        for start, end in zip(self.control_points, self.control_points[1:]):
            delta = (
                end.cumulative_meterset_weight
                - start.cumulative_meterset_weight
            )
            result.append(
                DeliverySegment(
                    start_control_point_index=start.index,
                    end_control_point_index=end.index,
                    delta_cumulative_meterset_weight=delta,
                    delta_mu=self.beam_meterset_mu * delta / final_weight,
                )
            )
        return tuple(result)

    def device_definition(self, device_type: str) -> BeamLimitingDeviceDefinition:
        for item in self.device_definitions:
            if item.device_type == device_type:
                return item
        raise KeyError(device_type)


@dataclass(frozen=True)
class RtPlan:
    sop_instance_uid: str
    plan_label: str | None
    plan_name: str | None
    fraction_group_number: int
    number_of_fractions_planned: int | None
    beams: tuple[Beam, ...]

    def beam(self, number: int) -> Beam:
        for item in self.beams:
            if item.number == number:
                return item
        raise KeyError(number)


def _optional_float(dataset: pydicom.dataset.Dataset, name: str) -> float | None:
    value = getattr(dataset, name, None)
    return None if value is None or value == "" else float(value)


def _optional_str(dataset: pydicom.dataset.Dataset, name: str) -> str | None:
    value = getattr(dataset, name, None)
    return None if value is None or str(value) == "" else str(value)


def _optional_float_tuple(
    dataset: pydicom.dataset.Dataset,
    name: str,
    length: int,
) -> tuple[float, ...] | None:
    value = getattr(dataset, name, None)
    if value is None:
        return None
    result = tuple(float(item) for item in value)
    if len(result) != length:
        raise ValueError(f"{name} must contain {length} values.")
    return result


def _inherit_float(
    dataset: pydicom.dataset.Dataset,
    name: str,
    previous: float | None,
) -> float | None:
    value = _optional_float(dataset, name)
    return previous if value is None else value


def _inherit_str(
    dataset: pydicom.dataset.Dataset,
    name: str,
    previous: str | None,
) -> str | None:
    value = _optional_str(dataset, name)
    return previous if value is None else value


def _inherit_tuple3(
    dataset: pydicom.dataset.Dataset,
    name: str,
    previous: tuple[float, float, float] | None,
) -> tuple[float, float, float] | None:
    value = _optional_float_tuple(dataset, name, 3)
    if value is None:
        return previous
    return (value[0], value[1], value[2])


def _select_fraction_group(
    dataset: pydicom.dataset.Dataset,
    requested_number: int | None,
) -> pydicom.dataset.Dataset:
    sequence = getattr(dataset, "FractionGroupSequence", None)
    if not sequence:
        raise ValueError("RTPLAN has no FractionGroupSequence.")

    if requested_number is None:
        if len(sequence) != 1:
            numbers = [int(item.FractionGroupNumber) for item in sequence]
            raise ValueError(
                "RTPLAN contains multiple fraction groups; select one explicitly: "
                + ", ".join(str(number) for number in numbers)
            )
        return sequence[0]

    for item in sequence:
        if int(item.FractionGroupNumber) == requested_number:
            return item
    raise ValueError(f"FractionGroupNumber {requested_number} not found.")


def _beam_metersets(
    fraction_group: pydicom.dataset.Dataset,
) -> dict[int, float]:
    referenced = getattr(fraction_group, "ReferencedBeamSequence", None)
    if not referenced:
        raise ValueError("Selected fraction group has no ReferencedBeamSequence.")

    result: dict[int, float] = {}
    for item in referenced:
        number = int(item.ReferencedBeamNumber)
        if number in result:
            raise ValueError(f"Beam {number} is referenced more than once.")
        if not hasattr(item, "BeamMeterset"):
            raise ValueError(f"Referenced beam {number} has no BeamMeterset.")
        result[number] = float(item.BeamMeterset)
    return result


def _device_definitions(
    beam: pydicom.dataset.Dataset,
) -> tuple[BeamLimitingDeviceDefinition, ...]:
    sequence = getattr(beam, "BeamLimitingDeviceSequence", None)
    if not sequence:
        raise ValueError(f"Beam {beam.BeamNumber} has no BeamLimitingDeviceSequence.")

    result: list[BeamLimitingDeviceDefinition] = []
    seen: set[str] = set()

    for item in sequence:
        device_type = str(item.RTBeamLimitingDeviceType)
        if device_type in seen:
            raise ValueError(
                f"Beam {beam.BeamNumber} defines device {device_type} more than once."
            )
        seen.add(device_type)

        pair_count = int(item.NumberOfLeafJawPairs)
        if pair_count < 1:
            raise ValueError(f"{device_type}: NumberOfLeafJawPairs must be positive.")

        boundaries = _optional_float_tuple(
            item,
            "LeafPositionBoundaries",
            pair_count + 1,
        )
        if device_type.startswith("MLC") and boundaries is None:
            raise ValueError(
                f"Beam {beam.BeamNumber} {device_type} has no LeafPositionBoundaries."
            )

        result.append(
            BeamLimitingDeviceDefinition(
                device_type=device_type,
                number_of_leaf_jaw_pairs=pair_count,
                leaf_position_boundaries_mm=boundaries,
                source_to_device_distance_mm=_optional_float(
                    item,
                    "SourceToBeamLimitingDeviceDistance",
                ),
            )
        )

    return tuple(result)


def _fluence_mode(
    beam: pydicom.dataset.Dataset,
) -> tuple[str | None, str | None]:
    sequence = getattr(beam, "PrimaryFluenceModeSequence", None)
    if not sequence:
        return None, None
    if len(sequence) != 1:
        raise ValueError(
            f"Beam {beam.BeamNumber} has multiple PrimaryFluenceModeSequence items."
        )
    item = sequence[0]
    return _optional_str(item, "FluenceMode"), _optional_str(item, "FluenceModeID")


def _control_points(
    beam: pydicom.dataset.Dataset,
    definitions: tuple[BeamLimitingDeviceDefinition, ...],
    *,
    cmw_tolerance: float,
) -> tuple[ControlPoint, ...]:
    sequence = getattr(beam, "ControlPointSequence", None)
    if not sequence:
        raise ValueError(f"Beam {beam.BeamNumber} has no ControlPointSequence.")

    expected_count = int(getattr(beam, "NumberOfControlPoints", len(sequence)))
    if expected_count != len(sequence):
        raise ValueError(
            f"Beam {beam.BeamNumber}: NumberOfControlPoints={expected_count}, "
            f"but sequence contains {len(sequence)} items."
        )

    definition_by_type = {item.device_type: item for item in definitions}
    current_devices: dict[str, tuple[float, ...]] = {}

    gantry: float | None = None
    gantry_direction: str | None = None
    collimator: float | None = None
    collimator_direction: str | None = None
    support: float | None = None
    support_direction: str | None = None
    eccentric: float | None = None
    eccentric_direction: str | None = None
    table_vertical: float | None = None
    table_longitudinal: float | None = None
    table_lateral: float | None = None
    energy: float | None = None
    dose_rate: float | None = None
    ssd: float | None = None
    isocenter: tuple[float, float, float] | None = None

    result: list[ControlPoint] = []

    for sequence_index, cp in enumerate(sequence):
        cp_index = int(cp.ControlPointIndex)
        if cp_index != sequence_index:
            raise ValueError(
                f"Beam {beam.BeamNumber}: expected ControlPointIndex "
                f"{sequence_index}, got {cp_index}."
            )

        cmw = float(cp.CumulativeMetersetWeight)

        positions_sequence = getattr(cp, "BeamLimitingDevicePositionSequence", None)
        if positions_sequence:
            for position_item in positions_sequence:
                device_type = str(position_item.RTBeamLimitingDeviceType)
                if device_type not in definition_by_type:
                    raise ValueError(
                        f"Beam {beam.BeamNumber} CP {cp_index} contains undefined "
                        f"device {device_type}."
                    )
                definition = definition_by_type[device_type]
                values = tuple(float(value) for value in position_item.LeafJawPositions)
                expected = 2 * definition.number_of_leaf_jaw_pairs
                if len(values) != expected:
                    raise ValueError(
                        f"Beam {beam.BeamNumber} CP {cp_index} {device_type}: "
                        f"expected {expected} positions, got {len(values)}."
                    )
                current_devices[device_type] = values

        missing_devices = [
            item.device_type
            for item in definitions
            if item.device_type not in current_devices
        ]
        if missing_devices:
            raise ValueError(
                f"Beam {beam.BeamNumber} CP {cp_index} has no inherited state for "
                + ", ".join(missing_devices)
            )

        gantry = _inherit_float(cp, "GantryAngle", gantry)
        gantry_direction = _inherit_str(
            cp, "GantryRotationDirection", gantry_direction
        )
        collimator = _inherit_float(cp, "BeamLimitingDeviceAngle", collimator)
        collimator_direction = _inherit_str(
            cp,
            "BeamLimitingDeviceRotationDirection",
            collimator_direction,
        )
        support = _inherit_float(cp, "PatientSupportAngle", support)
        support_direction = _inherit_str(
            cp, "PatientSupportRotationDirection", support_direction
        )
        eccentric = _inherit_float(cp, "TableTopEccentricAngle", eccentric)
        eccentric_direction = _inherit_str(
            cp,
            "TableTopEccentricRotationDirection",
            eccentric_direction,
        )
        table_vertical = _inherit_float(
            cp, "TableTopVerticalPosition", table_vertical
        )
        table_longitudinal = _inherit_float(
            cp, "TableTopLongitudinalPosition", table_longitudinal
        )
        table_lateral = _inherit_float(
            cp, "TableTopLateralPosition", table_lateral
        )
        energy = _inherit_float(cp, "NominalBeamEnergy", energy)
        dose_rate = _inherit_float(cp, "DoseRateSet", dose_rate)
        ssd = _inherit_float(cp, "SourceToSurfaceDistance", ssd)
        isocenter = _inherit_tuple3(cp, "IsocenterPosition", isocenter)

        states = tuple(
            BeamLimitingDeviceState(
                device_type=item.device_type,
                positions_mm=current_devices[item.device_type],
            )
            for item in definitions
        )

        result.append(
            ControlPoint(
                index=cp_index,
                cumulative_meterset_weight=cmw,
                gantry_angle_deg=gantry,
                gantry_rotation_direction=gantry_direction,
                collimator_angle_deg=collimator,
                collimator_rotation_direction=collimator_direction,
                patient_support_angle_deg=support,
                patient_support_rotation_direction=support_direction,
                table_top_eccentric_angle_deg=eccentric,
                table_top_eccentric_rotation_direction=eccentric_direction,
                table_top_vertical_position_mm=table_vertical,
                table_top_longitudinal_position_mm=table_longitudinal,
                table_top_lateral_position_mm=table_lateral,
                nominal_beam_energy_mv=energy,
                dose_rate_set_mu_min=dose_rate,
                source_to_surface_distance_mm=ssd,
                isocenter_position_mm=isocenter,
                device_positions=states,
            )
        )

    weights = np.asarray(
        [item.cumulative_meterset_weight for item in result],
        dtype=np.float64,
    )
    if abs(weights[0]) > cmw_tolerance:
        raise ValueError(
            f"Beam {beam.BeamNumber}: first CumulativeMetersetWeight is "
            f"{weights[0]}, expected 0."
        )
    if np.any(np.diff(weights) < -cmw_tolerance):
        raise ValueError(
            f"Beam {beam.BeamNumber}: CumulativeMetersetWeight is not monotonic."
        )

    return tuple(result)


def _parse_beam(
    beam: pydicom.dataset.Dataset,
    *,
    beam_meterset_mu: float,
    cmw_tolerance: float,
) -> Beam:
    definitions = _device_definitions(beam)
    control_points = _control_points(
        beam,
        definitions,
        cmw_tolerance=cmw_tolerance,
    )

    final_weight = _optional_float(beam, "FinalCumulativeMetersetWeight")
    if final_weight is None:
        final_weight = control_points[-1].cumulative_meterset_weight
    if final_weight <= 0:
        raise ValueError(
            f"Beam {beam.BeamNumber} has invalid FinalCumulativeMetersetWeight "
            f"{final_weight}."
        )

    last_weight = control_points[-1].cumulative_meterset_weight
    if not np.isclose(
        last_weight,
        final_weight,
        rtol=0.0,
        atol=cmw_tolerance,
    ):
        raise ValueError(
            f"Beam {beam.BeamNumber}: last CumulativeMetersetWeight "
            f"{last_weight} != FinalCumulativeMetersetWeight {final_weight}."
        )

    fluence_mode, fluence_mode_id = _fluence_mode(beam)

    return Beam(
        number=int(beam.BeamNumber),
        name=str(getattr(beam, "BeamName", f"Beam{beam.BeamNumber}")),
        description=_optional_str(beam, "BeamDescription"),
        beam_type=str(beam.BeamType),
        radiation_type=str(beam.RadiationType),
        treatment_delivery_type=str(beam.TreatmentDeliveryType),
        treatment_machine_name=_optional_str(beam, "TreatmentMachineName"),
        source_axis_distance_mm=_optional_float(beam, "SourceAxisDistance"),
        final_cumulative_meterset_weight=float(final_weight),
        beam_meterset_mu=float(beam_meterset_mu),
        fluence_mode=fluence_mode,
        fluence_mode_id=fluence_mode_id,
        device_definitions=definitions,
        control_points=control_points,
    )


def load_rtplan(
    path: str | Path,
    *,
    fraction_group_number: int | None = None,
    cmw_tolerance: float = 1e-6,
) -> RtPlan:
    """Load one photon RTPLAN into an explicit beam/control-point model.

    Patient identity is deliberately not retained. BeamMeterset is read from
    the selected FractionGroupSequence and matched by ReferencedBeamNumber.

    Control-point attributes and beam-limiting-device positions follow DICOM
    state inheritance: values omitted at a later control point retain the
    previous explicitly supplied state.
    """

    dataset = pydicom.dcmread(path, force=False)
    if getattr(dataset, "Modality", None) != "RTPLAN":
        raise ValueError("Input file is not an RT Plan.")

    beam_sequence = getattr(dataset, "BeamSequence", None)
    if not beam_sequence:
        raise ValueError("RTPLAN has no BeamSequence.")

    fraction_group = _select_fraction_group(dataset, fraction_group_number)
    metersets = _beam_metersets(fraction_group)

    beam_by_number = {int(item.BeamNumber): item for item in beam_sequence}
    if len(beam_by_number) != len(beam_sequence):
        raise ValueError("BeamNumber values are not unique.")

    parsed: list[Beam] = []
    for number, meterset in metersets.items():
        if number not in beam_by_number:
            raise ValueError(
                f"Fraction group references BeamNumber {number}, "
                "but BeamSequence does not contain it."
            )
        parsed.append(
            _parse_beam(
                beam_by_number[number],
                beam_meterset_mu=meterset,
                cmw_tolerance=cmw_tolerance,
            )
        )

    return RtPlan(
        sop_instance_uid=str(dataset.SOPInstanceUID),
        plan_label=_optional_str(dataset, "RTPlanLabel"),
        plan_name=_optional_str(dataset, "RTPlanName"),
        fraction_group_number=int(fraction_group.FractionGroupNumber),
        number_of_fractions_planned=(
            int(fraction_group.NumberOfFractionsPlanned)
            if hasattr(fraction_group, "NumberOfFractionsPlanned")
            else None
        ),
        beams=tuple(parsed),
    )
