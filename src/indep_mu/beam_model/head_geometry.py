from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml


@dataclass(frozen=True)
class HeadGeometryMessage:
    severity: str
    code: str
    message: str


@dataclass(frozen=True)
class RoundedLeafTipGeometry:
    """Physical parameters required by BEAMnrc SYNCMLCE ENDTYPE=0."""

    zmin_cm: float | None
    zmax_cm: float | None
    leaf_radius_cm: float | None
    cylinder_axis_z_cm: float | None
    leaf_bank_tilt_rad: float | None
    negative_bank: int | None


@dataclass(frozen=True)
class FocusedJawGeometry:
    device_type: str
    zmin_cm: float | None
    zmax_cm: float | None
    negative_bank: int | None


@dataclass(frozen=True)
class VersaHdHeadGeometry:
    machine_model_id: str
    mlc_end_type: str
    mlc: RoundedLeafTipGeometry
    jaws: tuple[FocusedJawGeometry, ...]
    evidence: tuple[str, ...] = ()


@dataclass(frozen=True)
class HeadGeometryReadiness:
    messages: tuple[HeadGeometryMessage, ...]

    @property
    def errors(self) -> tuple[HeadGeometryMessage, ...]:
        return tuple(item for item in self.messages if item.severity == "ERROR")

    @property
    def warnings(self) -> tuple[HeadGeometryMessage, ...]:
        return tuple(item for item in self.messages if item.severity == "WARNING")

    @property
    def ready_for_source_focused_jaws(self) -> bool:
        return not any(
            item.severity == "ERROR" and item.code.startswith("JAW_")
            for item in self.messages
        )

    @property
    def ready_for_rounded_agility_syncmlce(self) -> bool:
        # Geometry completeness alone is not enough. The DICOM projected
        # leaf-edge -> BEAMnrc cylinder-origin transform still needs a
        # validated implementation.
        return False


def _optional_float(raw: dict[str, Any], key: str) -> float | None:
    value = raw.get(key)
    return None if value is None else float(value)


def _optional_int(raw: dict[str, Any], key: str) -> int | None:
    value = raw.get(key)
    return None if value is None else int(value)


def load_head_geometry_yaml(path: str | Path) -> VersaHdHeadGeometry:
    raw = yaml.safe_load(Path(path).read_text(encoding="utf-8"))
    if not isinstance(raw, dict):
        raise ValueError("Head geometry YAML root must be a mapping.")

    mlc_raw = raw.get("mlc")
    if not isinstance(mlc_raw, dict):
        raise ValueError("mlc mapping is required.")

    jaws_raw = raw.get("jaws", [])
    if not isinstance(jaws_raw, list):
        raise ValueError("jaws must be a list.")

    evidence_raw = raw.get("evidence", [])
    if not isinstance(evidence_raw, list):
        raise ValueError("evidence must be a list.")

    return VersaHdHeadGeometry(
        machine_model_id=str(raw.get("machine_model_id", "")).strip(),
        mlc_end_type=str(mlc_raw.get("end_type", "")).strip().lower(),
        mlc=RoundedLeafTipGeometry(
            zmin_cm=_optional_float(mlc_raw, "zmin_cm"),
            zmax_cm=_optional_float(mlc_raw, "zmax_cm"),
            leaf_radius_cm=_optional_float(mlc_raw, "leaf_radius_cm"),
            cylinder_axis_z_cm=_optional_float(mlc_raw, "cylinder_axis_z_cm"),
            leaf_bank_tilt_rad=_optional_float(mlc_raw, "leaf_bank_tilt_rad"),
            negative_bank=_optional_int(mlc_raw, "negative_bank"),
        ),
        jaws=tuple(
            FocusedJawGeometry(
                device_type=str(item.get("device_type", "")).strip(),
                zmin_cm=_optional_float(item, "zmin_cm"),
                zmax_cm=_optional_float(item, "zmax_cm"),
                negative_bank=_optional_int(item, "negative_bank"),
            )
            for item in jaws_raw
            if isinstance(item, dict)
        ),
        evidence=tuple(str(item).strip() for item in evidence_raw if str(item).strip()),
    )


def validate_head_geometry(
    geometry: VersaHdHeadGeometry,
) -> HeadGeometryReadiness:
    """Check whether physical head inputs are sufficiently specified.

    This validator intentionally distinguishes "numbers are present" from
    "the rounded Agility mapping is physically validated". Even a complete
    parameter set cannot enable production ENDTYPE=0 mapping until the
    cylinder-origin transform has been independently verified.
    """

    messages: list[HeadGeometryMessage] = []

    if not geometry.machine_model_id:
        messages.append(
            HeadGeometryMessage("ERROR", "MACHINE_ID", "machine_model_id is missing.")
        )

    if geometry.mlc_end_type != "rounded":
        messages.append(
            HeadGeometryMessage(
                "ERROR",
                "MLC_END_TYPE",
                "Agility production model is expected to use an explicitly "
                "validated rounded-tip representation.",
            )
        )

    required_mlc = {
        "zmin_cm": geometry.mlc.zmin_cm,
        "zmax_cm": geometry.mlc.zmax_cm,
        "leaf_radius_cm": geometry.mlc.leaf_radius_cm,
        "cylinder_axis_z_cm": geometry.mlc.cylinder_axis_z_cm,
        "leaf_bank_tilt_rad": geometry.mlc.leaf_bank_tilt_rad,
        "negative_bank": geometry.mlc.negative_bank,
    }
    for name, value in required_mlc.items():
        if value is None:
            messages.append(
                HeadGeometryMessage(
                    "ERROR",
                    f"MLC_{name.upper()}",
                    f"Physical Agility parameter {name} is not supplied.",
                )
            )

    if (
        geometry.mlc.zmin_cm is not None
        and geometry.mlc.zmax_cm is not None
        and not (0 < geometry.mlc.zmin_cm < geometry.mlc.zmax_cm)
    ):
        messages.append(
            HeadGeometryMessage(
                "ERROR",
                "MLC_Z_RANGE",
                "MLC geometry requires 0 < zmin_cm < zmax_cm.",
            )
        )
    if geometry.mlc.leaf_radius_cm is not None and geometry.mlc.leaf_radius_cm <= 0:
        messages.append(
            HeadGeometryMessage(
                "ERROR",
                "MLC_RADIUS",
                "leaf_radius_cm must be positive.",
            )
        )
    if geometry.mlc.negative_bank is not None and geometry.mlc.negative_bank not in {1, 2}:
        messages.append(
            HeadGeometryMessage(
                "ERROR",
                "MLC_NEGATIVE_BANK",
                "negative_bank must be 1 or 2.",
            )
        )

    if not geometry.jaws:
        messages.append(
            HeadGeometryMessage(
                "ERROR",
                "JAW_MISSING",
                "No physical jaw geometry has been supplied.",
            )
        )

    seen_devices: set[str] = set()
    for jaw in geometry.jaws:
        prefix = f"JAW_{jaw.device_type or 'UNKNOWN'}"
        if not jaw.device_type:
            messages.append(
                HeadGeometryMessage("ERROR", f"{prefix}_TYPE", "Jaw device_type is empty.")
            )
        elif jaw.device_type in seen_devices:
            messages.append(
                HeadGeometryMessage(
                    "ERROR",
                    f"{prefix}_DUPLICATE",
                    f"Jaw device {jaw.device_type} is defined more than once.",
                )
            )
        seen_devices.add(jaw.device_type)

        if jaw.zmin_cm is None or jaw.zmax_cm is None:
            messages.append(
                HeadGeometryMessage(
                    "ERROR",
                    f"{prefix}_Z",
                    f"{jaw.device_type or 'Jaw'} requires zmin_cm and zmax_cm.",
                )
            )
        elif not (0 < jaw.zmin_cm < jaw.zmax_cm):
            messages.append(
                HeadGeometryMessage(
                    "ERROR",
                    f"{prefix}_Z",
                    f"{jaw.device_type} requires 0 < zmin_cm < zmax_cm.",
                )
            )

        if jaw.negative_bank not in {1, 2}:
            messages.append(
                HeadGeometryMessage(
                    "ERROR",
                    f"{prefix}_NEGATIVE_BANK",
                    f"{jaw.device_type or 'Jaw'} negative_bank must be 1 or 2.",
                )
            )

    messages.append(
        HeadGeometryMessage(
            "ERROR",
            "MLC_ROUNDED_MAPPING_UNVALIDATED",
            "DICOM projected Agility leaf position -> SYNCMLCE ENDTYPE=0 "
            "cylinder-origin transform is not yet independently validated.",
        )
    )

    if not geometry.evidence:
        messages.append(
            HeadGeometryMessage(
                "WARNING",
                "GEOMETRY_EVIDENCE",
                "No provenance/evidence entries are attached to the head geometry.",
            )
        )

    return HeadGeometryReadiness(messages=tuple(messages))
