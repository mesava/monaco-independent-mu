import pytest

from indep_mu.patient_model.compositions import BASE_COMPOSITIONS, mix_mass_fractions


@pytest.mark.parametrize("name", list(BASE_COMPOSITIONS))
def test_reference_compositions_are_normalized(name: str) -> None:
    composition = BASE_COMPOSITIONS[name]
    assert sum(composition.mass_fractions.values()) == pytest.approx(1.0, abs=1e-6)


def test_mixture_is_normalized() -> None:
    mixed = mix_mass_fractions(
        BASE_COMPOSITIONS["MuscleSkeletalIcrp"],
        BASE_COMPOSITIONS["BoneCorticalIcrp"],
        0.75,
    )
    assert sum(mixed.values()) == pytest.approx(1.0, abs=1e-12)
