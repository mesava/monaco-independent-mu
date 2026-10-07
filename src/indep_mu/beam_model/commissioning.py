from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import yaml


@dataclass(frozen=True)
class FieldSize:
    x_mm: float
    y_mm: float

    def __post_init__(self) -> None:
        if self.x_mm <= 0 or self.y_mm <= 0:
            raise ValueError("Field dimensions must be positive.")

    @property
    def equivalent_square_mm(self) -> float:
        return 2.0 * self.x_mm * self.y_mm / (self.x_mm + self.y_mm)


@dataclass(frozen=True)
class ReferenceCalibration:
    """Measured absolute calibration with explicit delivered MU.

    The explicit delivered_mu field is mandatory by design. A value such as
    0.995 Gy is never interpreted as Gy/MU or Gy/100 MU without this metadata.
    """

    field: FieldSize
    ssd_mm: float
    depth_mm: float
    delivered_mu: float
    measured_dose_gy: float

    def __post_init__(self) -> None:
        if self.ssd_mm <= 0:
            raise ValueError("SSD must be positive.")
        if self.depth_mm < 0:
            raise ValueError("Depth cannot be negative.")
        if self.delivered_mu <= 0:
            raise ValueError("delivered_mu must be positive.")
        if self.measured_dose_gy <= 0:
            raise ValueError("measured_dose_gy must be positive.")

    @property
    def dose_gy_per_mu(self) -> float:
        return self.measured_dose_gy / self.delivered_mu

    @property
    def dose_gy_per_100mu(self) -> float:
        return 100.0 * self.dose_gy_per_mu


@dataclass(frozen=True)
class OutputFactorPoint:
    field: FieldSize
    output_factor: float

    def __post_init__(self) -> None:
        if self.output_factor <= 0:
            raise ValueError("Output factor must be positive.")


@dataclass(frozen=True)
class DepthDoseCurve:
    field: FieldSize
    ssd_mm: float
    depth_mm: np.ndarray
    relative_dose: np.ndarray
    normalization: str = "DMAX"

    def __post_init__(self) -> None:
        depth = np.asarray(self.depth_mm, dtype=np.float64)
        dose = np.asarray(self.relative_dose, dtype=np.float64)
        if depth.ndim != 1 or dose.ndim != 1 or depth.size != dose.size:
            raise ValueError("PDD depth and dose arrays must be equal 1-D arrays.")
        if depth.size < 2:
            raise ValueError("PDD curve requires at least two samples.")
        if np.any(np.diff(depth) <= 0):
            raise ValueError("PDD depths must be strictly increasing.")
        if np.any(dose < 0) or not np.all(np.isfinite(dose)):
            raise ValueError("PDD values must be finite and non-negative.")
        if self.ssd_mm <= 0:
            raise ValueError("PDD SSD must be positive.")


@dataclass(frozen=True)
class BeamProfile:
    field: FieldSize
    ssd_mm: float
    depth_mm: float
    axis: str
    position_mm: np.ndarray
    relative_dose: np.ndarray

    def __post_init__(self) -> None:
        axis = self.axis.upper()
        if axis not in {"INPLANE", "CROSSPLANE"}:
            raise ValueError("Profile axis must be INPLANE or CROSSPLANE.")
        position = np.asarray(self.position_mm, dtype=np.float64)
        dose = np.asarray(self.relative_dose, dtype=np.float64)
        if position.ndim != 1 or dose.ndim != 1 or position.size != dose.size:
            raise ValueError("Profile position and dose must be equal 1-D arrays.")
        if position.size < 2:
            raise ValueError("Profile requires at least two samples.")
        if np.any(np.diff(position) <= 0):
            raise ValueError("Profile positions must be strictly increasing.")
        if np.any(dose < 0) or not np.all(np.isfinite(dose)):
            raise ValueError("Profile dose must be finite and non-negative.")
        if self.ssd_mm <= 0 or self.depth_mm < 0:
            raise ValueError("Invalid profile geometry.")


@dataclass(frozen=True)
class EnergyCommissioningData:
    energy_model_id: str
    calibration: ReferenceCalibration
    output_factors: tuple[OutputFactorPoint, ...] = ()
    pdd_curves: tuple[DepthDoseCurve, ...] = ()
    profiles: tuple[BeamProfile, ...] = ()

    def validate_reference_output_factor(
        self,
        *,
        reference_field_mm: float = 100.0,
        tolerance: float = 0.02,
    ) -> None:
        matching = [
            point
            for point in self.output_factors
            if abs(point.field.x_mm - reference_field_mm) < 1e-9
            and abs(point.field.y_mm - reference_field_mm) < 1e-9
        ]
        if not matching:
            raise ValueError(
                f"No {reference_field_mm:g} x {reference_field_mm:g} mm "
                "reference output factor is present."
            )
        if len(matching) > 1:
            raise ValueError("Reference output factor is duplicated.")
        value = matching[0].output_factor
        if abs(value - 1.0) > tolerance:
            raise ValueError(
                f"Reference output factor {value} differs from 1 by more "
                f"than tolerance {tolerance}."
            )


@dataclass(frozen=True)
class CommissioningDataset:
    machine_model_id: str
    energies: tuple[EnergyCommissioningData, ...]

    def __post_init__(self) -> None:
        ids = [item.energy_model_id for item in self.energies]
        if len(ids) != len(set(ids)):
            raise ValueError("Commissioning energy_model_id values must be unique.")

    def energy(self, model_id: str) -> EnergyCommissioningData:
        for item in self.energies:
            if item.energy_model_id == model_id:
                return item
        raise KeyError(model_id)


def _field(data: dict[str, Any]) -> FieldSize:
    try:
        return FieldSize(
            x_mm=float(data["x_mm"]),
            y_mm=float(data["y_mm"]),
        )
    except KeyError as exc:
        raise ValueError(f"Field definition is missing {exc.args[0]!r}.") from exc


def _curve_array(data: dict[str, Any], key: str) -> np.ndarray:
    if key not in data:
        raise ValueError(f"Curve definition is missing {key!r}.")
    return np.asarray([float(value) for value in data[key]], dtype=np.float64)


def load_commissioning_yaml(path: str | Path) -> CommissioningDataset:
    """Load a version-controlled commissioning dataset from YAML."""

    source = Path(path)
    raw = yaml.safe_load(source.read_text(encoding="utf-8"))
    if not isinstance(raw, dict):
        raise ValueError("Commissioning YAML root must be a mapping.")

    machine_model_id = str(raw.get("machine_model_id", "")).strip()
    if not machine_model_id:
        raise ValueError("machine_model_id is required.")

    raw_energies = raw.get("energies")
    if not isinstance(raw_energies, list) or not raw_energies:
        raise ValueError("energies must be a non-empty list.")

    energies: list[EnergyCommissioningData] = []
    for item in raw_energies:
        if not isinstance(item, dict):
            raise ValueError("Each energy entry must be a mapping.")

        model_id = str(item.get("energy_model_id", "")).strip()
        if not model_id:
            raise ValueError("energy_model_id is required.")

        calibration_data = item.get("calibration")
        if not isinstance(calibration_data, dict):
            raise ValueError(f"{model_id}: calibration mapping is required.")
        if "delivered_mu" not in calibration_data:
            raise ValueError(
                f"{model_id}: calibration.delivered_mu is mandatory; "
                "the loader never assumes 100 MU."
            )

        calibration = ReferenceCalibration(
            field=_field(calibration_data["field"]),
            ssd_mm=float(calibration_data["ssd_mm"]),
            depth_mm=float(calibration_data["depth_mm"]),
            delivered_mu=float(calibration_data["delivered_mu"]),
            measured_dose_gy=float(calibration_data["measured_dose_gy"]),
        )

        output_factors = tuple(
            OutputFactorPoint(
                field=_field(point["field"]),
                output_factor=float(point["output_factor"]),
            )
            for point in item.get("output_factors", [])
        )

        pdd_curves = tuple(
            DepthDoseCurve(
                field=_field(curve["field"]),
                ssd_mm=float(curve["ssd_mm"]),
                depth_mm=_curve_array(curve, "depth_mm"),
                relative_dose=_curve_array(curve, "relative_dose"),
                normalization=str(curve.get("normalization", "DMAX")),
            )
            for curve in item.get("pdd_curves", [])
        )

        profiles = tuple(
            BeamProfile(
                field=_field(profile["field"]),
                ssd_mm=float(profile["ssd_mm"]),
                depth_mm=float(profile["depth_mm"]),
                axis=str(profile["axis"]),
                position_mm=_curve_array(profile, "position_mm"),
                relative_dose=_curve_array(profile, "relative_dose"),
            )
            for profile in item.get("profiles", [])
        )

        energies.append(
            EnergyCommissioningData(
                energy_model_id=model_id,
                calibration=calibration,
                output_factors=output_factors,
                pdd_curves=pdd_curves,
                profiles=profiles,
            )
        )

    return CommissioningDataset(
        machine_model_id=machine_model_id,
        energies=tuple(energies),
    )
