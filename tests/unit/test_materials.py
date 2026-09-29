import pytest

from indep_mu.patient_model.materials import material_blend_for_density


def test_pure_air_interval() -> None:
    blend = material_blend_for_density(0.001)
    assert blend.lower_material == "DryAir"
    assert blend.upper_material == "DryAir"
    assert blend.is_pure


def test_muscle_to_bone_interpolation() -> None:
    blend = material_blend_for_density((1.08 + 1.85) / 2)
    assert blend.lower_material == "MuscleSkeletalIcrp"
    assert blend.upper_material == "BoneCorticalIcrp"
    assert blend.upper_fraction == pytest.approx(0.5)


def test_bone_interval_is_pure() -> None:
    blend = material_blend_for_density(2.57)
    assert blend.lower_material == "BoneCorticalIcrp"
    assert blend.upper_material == "BoneCorticalIcrp"
    assert blend.is_pure


def test_equal_threshold_uses_right_hand_material() -> None:
    blend = material_blend_for_density(3.0)
    assert blend.lower_material == "Titanium"
    assert blend.upper_material == "Titanium"
