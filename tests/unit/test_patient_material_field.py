import numpy as np

from indep_mu.patient_model.discretize import discretize_material_blends
from indep_mu.patient_model.materials import material_blend_arrays_for_density


def test_vectorized_material_assignment_across_patient_lut() -> None:
    rho = np.asarray([0.001, 0.5, 0.865, 0.93, 1.00, 1.465, 2.57, 3.0])
    field = material_blend_arrays_for_density(rho)

    names = field.material_names

    assert names[int(field.lower_material_id[0])] == "DryAir"
    assert names[int(field.upper_material_id[0])] == "DryAir"

    assert names[int(field.lower_material_id[1])] == "DryAir"
    assert names[int(field.upper_material_id[1])] == "MuscleSkeletalIcrp"

    assert names[int(field.lower_material_id[2])] == "MuscleSkeletalIcrp"
    assert names[int(field.upper_material_id[2])] == "AdiposeTissueICRP"

    assert names[int(field.lower_material_id[6])] == "BoneCorticalIcrp"
    assert names[int(field.upper_material_id[6])] == "BoneCorticalIcrp"

    # 3.0 g/cm3 is the abrupt Bone -> Titanium boundary: exact tie is Titanium.
    assert names[int(field.lower_material_id[7])] == "Titanium"
    assert names[int(field.upper_material_id[7])] == "Titanium"


def test_discretizer_canonicalizes_reverse_material_pair() -> None:
    # 0.865 lies halfway Muscle->Adipose, 0.995 halfway Adipose->Muscle.
    rho = np.asarray([0.865, 0.995])
    field = material_blend_arrays_for_density(rho)
    discrete = discretize_material_blends(field, mixture_bins=10)

    assert len(discrete.media) == 1
    medium = discrete.media[0]
    assert {medium.material_a, medium.material_b} == {
        "MuscleSkeletalIcrp",
        "AdiposeTissueICRP",
    }
    assert medium.fraction_b == 0.5


def test_egsphant_media_limit_is_enforced() -> None:
    rho = np.linspace(0.003, 2.99, 5000)
    field = material_blend_arrays_for_density(rho)

    # Artificially tiny cap to exercise the safety check.
    try:
        discretize_material_blends(field, mixture_bins=19, max_media=2)
    except ValueError as exc:
        assert "exceeding" in str(exc)
    else:
        raise AssertionError("Expected media-limit failure.")
