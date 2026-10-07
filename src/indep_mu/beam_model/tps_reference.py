from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import yaml

from .commissioning import FieldSize


@dataclass(frozen=True)
class TpsReferenceGeometry:
    field: FieldSize
    ssd_mm: float
    depth_mm: float
    delivered_mu: float

    def __post_init__(self) -> None:
        if self.ssd_mm <= 0:
            raise ValueError("SSD must be positive.")
        if self.depth_mm < 0:
            raise ValueError("Depth cannot be negative.")
        if self.delivered_mu <= 0:
            raise ValueError("delivered_mu must be positive.")


@dataclass(frozen=True)
class TpsReferenceDose:
    energy_model_id: str
    dose_gy: float

    def __post_init__(self) -> None:
        if not self.energy_model_id:
            raise ValueError("energy_model_id is required.")
        if self.dose_gy <= 0:
            raise ValueError("TPS reference dose must be positive.")


@dataclass(frozen=True)
class TpsReferenceDataset:
    machine_model_id: str
    tps_name: str
    tps_version: str | None
    geometry: TpsReferenceGeometry
    energies: tuple[TpsReferenceDose, ...]

    def __post_init__(self) -> None:
        ids = [item.energy_model_id for item in self.energies]
        if len(ids) != len(set(ids)):
            raise ValueError("TPS reference energy_model_id values must be unique.")

    def energy(self, model_id: str) -> TpsReferenceDose:
        for item in self.energies:
            if item.energy_model_id == model_id:
                return item
        raise KeyError(model_id)

    def dose_gy_per_mu(self, model_id: str) -> float:
        return self.energy(model_id).dose_gy / self.geometry.delivered_mu

    def dose_gy_per_100mu(self, model_id: str) -> float:
        return 100.0 * self.dose_gy_per_mu(model_id)


def load_tps_reference_yaml(path: str | Path) -> TpsReferenceDataset:
    raw = yaml.safe_load(Path(path).read_text(encoding="utf-8"))
    if not isinstance(raw, dict):
        raise ValueError("TPS reference YAML root must be a mapping.")

    geometry_raw = raw.get("reference_geometry")
    if not isinstance(geometry_raw, dict):
        raise ValueError("reference_geometry mapping is required.")
    field_raw = geometry_raw.get("field")
    if not isinstance(field_raw, dict):
        raise ValueError("reference_geometry.field mapping is required.")
    if "delivered_mu" not in geometry_raw:
        raise ValueError("reference_geometry.delivered_mu is mandatory.")

    energies_raw = raw.get("energies")
    if not isinstance(energies_raw, list) or not energies_raw:
        raise ValueError("energies must be a non-empty list.")

    dataset = TpsReferenceDataset(
        machine_model_id=str(raw.get("machine_model_id", "")).strip(),
        tps_name=str(raw.get("tps_name", "")).strip(),
        tps_version=(
            str(raw["tps_version"]).strip()
            if raw.get("tps_version") is not None
            else None
        ),
        geometry=TpsReferenceGeometry(
            field=FieldSize(
                x_mm=float(field_raw["x_mm"]),
                y_mm=float(field_raw["y_mm"]),
            ),
            ssd_mm=float(geometry_raw["ssd_mm"]),
            depth_mm=float(geometry_raw["depth_mm"]),
            delivered_mu=float(geometry_raw["delivered_mu"]),
        ),
        energies=tuple(
            TpsReferenceDose(
                energy_model_id=str(item["energy_model_id"]).strip(),
                dose_gy=float(item["dose_gy"]),
            )
            for item in energies_raw
        ),
    )

    if not dataset.machine_model_id:
        raise ValueError("machine_model_id is required.")
    if not dataset.tps_name:
        raise ValueError("tps_name is required.")

    return dataset
