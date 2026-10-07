from __future__ import annotations

from dataclasses import dataclass

from indep_mu.dicom.rtplan import Beam


@dataclass(frozen=True)
class EnergyModel:
    model_id: str
    nominal_energy_mv: float
    filter_mode: str  # "FF" or "FFF"
    dicom_fluence_mode: str | None = None
    dicom_fluence_mode_ids: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if self.filter_mode not in {"FF", "FFF"}:
            raise ValueError("filter_mode must be 'FF' or 'FFF'.")


@dataclass(frozen=True)
class MachineModel:
    machine_id: str
    manufacturer: str
    family: str
    mlc_family: str
    source_axis_distance_mm: float
    energy_models: tuple[EnergyModel, ...]

    def energy_model(self, model_id: str) -> EnergyModel:
        for item in self.energy_models:
            if item.model_id == model_id:
                return item
        raise KeyError(model_id)


VERSA_HD = MachineModel(
    machine_id="VERSA_HD",
    manufacturer="Elekta",
    family="Versa HD",
    mlc_family="Agility",
    source_axis_distance_mm=1000.0,
    energy_models=(
        EnergyModel("6MV", 6.0, "FF", dicom_fluence_mode="STANDARD"),
        EnergyModel(
            "6FFF",
            6.0,
            "FFF",
            dicom_fluence_mode="NON_STANDARD",
            dicom_fluence_mode_ids=("FFF", "6FFF", "6 FFF"),
        ),
        EnergyModel("10MV", 10.0, "FF", dicom_fluence_mode="STANDARD"),
        EnergyModel(
            "10FFF",
            10.0,
            "FFF",
            dicom_fluence_mode="NON_STANDARD",
            dicom_fluence_mode_ids=("FFF", "10FFF", "10 FFF"),
        ),
    ),
)


def _beam_energy(beam: Beam) -> float:
    energies = {
        float(cp.nominal_beam_energy_mv)
        for cp in beam.control_points
        if cp.nominal_beam_energy_mv is not None
    }
    if len(energies) != 1:
        raise ValueError(
            f"Beam {beam.number} must contain exactly one NominalBeamEnergy, "
            f"got {sorted(energies)}."
        )
    return next(iter(energies))


def resolve_energy_model(
    beam: Beam,
    machine: MachineModel = VERSA_HD,
    *,
    override_model_id: str | None = None,
    energy_tolerance_mv: float = 1e-6,
) -> EnergyModel:
    """Map DICOM beam metadata to one commissioned machine energy model.

    Energy alone is deliberately insufficient when both FF and FFF models exist
    at the same nominal MV. Ambiguous metadata is a hard error rather than an
    implicit guess.
    """

    energy = _beam_energy(beam)

    if override_model_id is not None:
        model = machine.energy_model(override_model_id)
        if abs(model.nominal_energy_mv - energy) > energy_tolerance_mv:
            raise ValueError(
                f"Override {override_model_id} nominal energy "
                f"{model.nominal_energy_mv} MV does not match DICOM {energy} MV."
            )
        return model

    candidates = [
        item
        for item in machine.energy_models
        if abs(item.nominal_energy_mv - energy) <= energy_tolerance_mv
    ]
    if not candidates:
        raise ValueError(
            f"No commissioned energy model for DICOM NominalBeamEnergy={energy} MV."
        )
    if len(candidates) == 1:
        return candidates[0]

    fluence_mode = beam.fluence_mode
    fluence_id = beam.fluence_mode_id

    if fluence_mode is not None:
        by_mode = [
            item
            for item in candidates
            if item.dicom_fluence_mode == fluence_mode
        ]
        if len(by_mode) == 1:
            return by_mode[0]
        if by_mode:
            candidates = by_mode

    if fluence_id is not None:
        key = fluence_id.casefold().replace("_", " ").strip()
        matched = []
        for item in candidates:
            aliases = {
                value.casefold().replace("_", " ").strip()
                for value in item.dicom_fluence_mode_ids
            }
            if key in aliases:
                matched.append(item)
        if len(matched) == 1:
            return matched[0]

    raise ValueError(
        f"Beam {beam.number}: DICOM metadata is insufficient to distinguish "
        f"commissioned models {[item.model_id for item in candidates]}. "
        "Provide an explicit mapping/override instead of inferring FF vs FFF."
    )
