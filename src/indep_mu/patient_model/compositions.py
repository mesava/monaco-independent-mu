from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class MaterialComposition:
    name: str
    reference_density_g_cm3: float
    mean_excitation_energy_ev: float
    mass_fractions: dict[str, float]
    source: str

    def __post_init__(self) -> None:
        total = sum(self.mass_fractions.values())
        if abs(total - 1.0) > 1e-6:
            raise ValueError(f"{self.name}: mass fractions sum to {total}, not 1.")


# Independent reference compositions. These are not extracted from Monaco.
# They follow NIST STAR entries with the same ICRP material names.
BASE_COMPOSITIONS: dict[str, MaterialComposition] = {
    "DryAir": MaterialComposition(
        name="DryAir",
        reference_density_g_cm3=1.205e-3,
        mean_excitation_energy_ev=85.7,
        mass_fractions={
            "C": 0.000124,
            "N": 0.755268,
            "O": 0.231781,
            "Ar": 0.012827,
        },
        source="NIST STAR: Air, Dry (near sea level)",
    ),
    "AdiposeTissueICRP": MaterialComposition(
        name="AdiposeTissueICRP",
        reference_density_g_cm3=0.920,
        mean_excitation_energy_ev=63.2,
        mass_fractions={
            "H": 0.119477,
            "C": 0.637240,
            "N": 0.007970,
            "O": 0.232333,
            "Na": 0.000500,
            "Mg": 0.000020,
            "P": 0.000160,
            "S": 0.000730,
            "Cl": 0.001190,
            "K": 0.000320,
            "Ca": 0.000020,
            "Fe": 0.000020,
            "Zn": 0.000020,
        },
        source="NIST STAR: Adipose Tissue (ICRP)",
    ),
    "MuscleSkeletalIcrp": MaterialComposition(
        name="MuscleSkeletalIcrp",
        reference_density_g_cm3=1.040,
        mean_excitation_energy_ev=75.3,
        mass_fractions={
            "H": 0.100637,
            "C": 0.107830,
            "N": 0.027680,
            "O": 0.754773,
            "Na": 0.000750,
            "Mg": 0.000190,
            "P": 0.001800,
            "S": 0.002410,
            "Cl": 0.000790,
            "K": 0.003020,
            "Ca": 0.000030,
            "Fe": 0.000040,
            "Zn": 0.000050,
        },
        source="NIST STAR: Muscle, Skeletal (ICRP)",
    ),
    "BoneCorticalIcrp": MaterialComposition(
        name="BoneCorticalIcrp",
        reference_density_g_cm3=1.850,
        mean_excitation_energy_ev=106.4,
        mass_fractions={
            "H": 0.047234,
            "C": 0.144330,
            "N": 0.041990,
            "O": 0.446096,
            "Mg": 0.002200,
            "P": 0.104970,
            "S": 0.003150,
            "Ca": 0.209930,
            "Zn": 0.000100,
        },
        source="NIST STAR: Bone, Cortical (ICRP)",
    ),
}


def mix_mass_fractions(
    material_a: MaterialComposition,
    material_b: MaterialComposition,
    fraction_b: float,
) -> dict[str, float]:
    """Linearly blend elemental mass fractions for a discrete mix medium."""

    if not 0.0 <= fraction_b <= 1.0:
        raise ValueError("fraction_b must be in [0, 1].")

    elements = set(material_a.mass_fractions) | set(material_b.mass_fractions)
    mixed = {
        element: (
            (1.0 - fraction_b) * material_a.mass_fractions.get(element, 0.0)
            + fraction_b * material_b.mass_fractions.get(element, 0.0)
        )
        for element in sorted(elements)
    }

    total = sum(mixed.values())
    return {element: value / total for element, value in mixed.items()}
