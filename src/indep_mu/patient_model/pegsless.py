from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Mapping

from .compositions import BASE_COMPOSITIONS, MaterialComposition, mix_mass_fractions
from .discretize import DiscreteMaterialField, DiscreteMedium


_ATOMIC_NUMBER = {
    "H": 1,
    "C": 6,
    "N": 7,
    "O": 8,
    "Na": 11,
    "Mg": 12,
    "P": 15,
    "S": 16,
    "Cl": 17,
    "Ar": 18,
    "K": 19,
    "Ca": 20,
    "Fe": 26,
    "Zn": 30,
}


@dataclass(frozen=True)
class PegslessEnergyBounds:
    """Optional cross-section energy bounds in MeV."""

    ae: float | None = None
    ue: float | None = None
    ap: float | None = None
    up: float | None = None

    def __post_init__(self) -> None:
        for name, value in (
            ("ae", self.ae),
            ("ue", self.ue),
            ("ap", self.ap),
            ("up", self.up),
        ):
            if value is not None and value <= 0:
                raise ValueError(f"{name} must be positive when supplied.")
        if self.ae is not None and self.ue is not None and self.ue <= self.ae:
            raise ValueError("ue must be greater than ae.")
        if self.ap is not None and self.up is not None and self.up <= self.ap:
            raise ValueError("up must be greater than ap.")


@dataclass(frozen=True)
class PegslessMedium:
    name: str
    elements: tuple[str, ...]
    mass_fractions: tuple[float, ...]
    bulk_density_g_cm3: float
    source_material_a: str
    source_material_b: str
    fraction_b: float

    def __post_init__(self) -> None:
        if not self.name or len(self.name) > 24:
            raise ValueError("PEGSless medium name must contain 1-24 characters.")
        if len(self.elements) != len(self.mass_fractions):
            raise ValueError("elements and mass_fractions must have equal length.")
        if not self.elements:
            raise ValueError("At least one element is required.")
        if self.bulk_density_g_cm3 <= 0:
            raise ValueError("bulk density must be positive.")
        if abs(sum(self.mass_fractions) - 1.0) > 1e-8:
            raise ValueError("PEGSless medium mass fractions must sum to 1.")
        if not 0.0 <= self.fraction_b <= 1.0:
            raise ValueError("fraction_b must be in [0, 1].")


@dataclass(frozen=True)
class PegslessMediaSet:
    media: tuple[PegslessMedium, ...]
    energy_bounds: PegslessEnergyBounds
    bremsstrahlung_correction: str = "NRC"

    def __post_init__(self) -> None:
        names = [medium.name for medium in self.media]
        if len(set(names)) != len(names):
            raise ValueError("PEGSless medium names must be unique.")
        if self.bremsstrahlung_correction not in {"KM", "NRC", "None"}:
            raise ValueError("Unsupported bremsstrahlung correction.")


def _stable_elements(mass_fractions: Mapping[str, float]) -> tuple[str, ...]:
    unknown = [element for element in mass_fractions if element not in _ATOMIC_NUMBER]
    if unknown:
        raise ValueError(
            "Atomic number is not defined for element(s): " + ", ".join(sorted(unknown))
        )
    return tuple(sorted(mass_fractions, key=lambda element: _ATOMIC_NUMBER[element]))


def _reference_density(
    material_a: MaterialComposition,
    material_b: MaterialComposition,
    fraction_b: float,
) -> float:
    """Reference density for cross-section generation.

    This is deliberately a composition-reference density, not the patient voxel
    density. DOSXYZnrc receives the actual voxel density from .egsphant and uses
    its RHOR/RHO(MEDIUM) scaling. The remaining density-effect sensitivity is
    quantified separately before commissioning is frozen.
    """

    return (
        (1.0 - fraction_b) * material_a.reference_density_g_cm3
        + fraction_b * material_b.reference_density_g_cm3
    )


def pegsless_medium_from_discrete(
    medium: DiscreteMedium,
    *,
    compositions: Mapping[str, MaterialComposition] = BASE_COMPOSITIONS,
) -> PegslessMedium:
    try:
        material_a = compositions[medium.material_a]
        material_b = compositions[medium.material_b]
    except KeyError as exc:
        raise ValueError(
            f"No independent elemental composition is defined for {exc.args[0]!r}."
        ) from exc

    if medium.material_a == medium.material_b or medium.fraction_b == 0.0:
        fractions = dict(material_a.mass_fractions)
        fraction_b = 0.0
    else:
        fractions = mix_mass_fractions(
            material_a,
            material_b,
            medium.fraction_b,
        )
        fraction_b = float(medium.fraction_b)

    elements = _stable_elements(fractions)
    values = tuple(float(fractions[element]) for element in elements)

    return PegslessMedium(
        name=medium.egsphant_name,
        elements=elements,
        mass_fractions=values,
        bulk_density_g_cm3=_reference_density(
            material_a,
            material_b,
            fraction_b,
        ),
        source_material_a=medium.material_a,
        source_material_b=medium.material_b,
        fraction_b=fraction_b,
    )


def build_pegsless_media_set(
    discrete: DiscreteMaterialField,
    *,
    energy_bounds: PegslessEnergyBounds | None = None,
    compositions: Mapping[str, MaterialComposition] = BASE_COMPOSITIONS,
    bremsstrahlung_correction: str = "NRC",
) -> PegslessMediaSet:
    """Build PEGSless definitions in exactly the .egsphant medium order."""

    media = tuple(
        pegsless_medium_from_discrete(item, compositions=compositions)
        for item in discrete.media
    )
    return PegslessMediaSet(
        media=media,
        energy_bounds=energy_bounds or PegslessEnergyBounds(),
        bremsstrahlung_correction=bremsstrahlung_correction,
    )


def render_media_definition(media_set: PegslessMediaSet) -> str:
    """Render an EGSnrc :start media definition: block.

    I-values are intentionally not injected for interpolated media. EGSnrc is
    allowed to calculate density-effect data on-the-fly from the elemental
    composition, as recommended for new pegsless materials in PIRS-701.
    """

    lines = [":start media definition:", ""]

    bounds = media_set.energy_bounds
    for key, value in (
        ("ae", bounds.ae),
        ("ue", bounds.ue),
        ("ap", bounds.ap),
        ("up", bounds.up),
    ):
        if value is not None:
            lines.append(f"    {key} = {value:.12g}")

    if any(value is not None for value in (bounds.ae, bounds.ue, bounds.ap, bounds.up)):
        lines.append("")

    for medium in media_set.media:
        lines.extend(
            [
                f"    :start {medium.name}:",
                f"        elements = {', '.join(medium.elements)}",
                "        mass fractions = "
                + ", ".join(f"{value:.12g}" for value in medium.mass_fractions),
                f"        bulk density = {medium.bulk_density_g_cm3:.12g}",
                "        bremsstrahlung correction = "
                + media_set.bremsstrahlung_correction,
                f"    :stop {medium.name}:",
                "",
            ]
        )

    lines.append(":stop media definition:")
    lines.append("")
    return "\n".join(lines)


def write_media_definition(path: str | Path, media_set: PegslessMediaSet) -> None:
    Path(path).write_text(render_media_definition(media_set), encoding="ascii")
