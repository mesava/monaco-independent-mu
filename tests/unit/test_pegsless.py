import numpy as np

from indep_mu.patient_model.discretize import discretize_material_blends
from indep_mu.patient_model.materials import material_blend_arrays_for_density
from indep_mu.patient_model.pegsless import (
    PegslessEnergyBounds,
    build_pegsless_media_set,
    render_media_definition,
)


def test_pegsless_media_follow_egsphant_order() -> None:
    density = np.asarray([0.001, 0.5, 0.92, 1.04, 1.5, 1.85])
    blends = material_blend_arrays_for_density(density)
    discrete = discretize_material_blends(blends, mixture_bins=5)

    media_set = build_pegsless_media_set(
        discrete,
        energy_bounds=PegslessEnergyBounds(
            ae=0.521,
            ue=20.511,
            ap=0.01,
            up=20.0,
        ),
    )

    assert [medium.name for medium in media_set.media] == [
        medium.egsphant_name for medium in discrete.media
    ]
    assert all(abs(sum(m.mass_fractions) - 1.0) < 1e-10 for m in media_set.media)


def test_render_media_definition_contains_required_pegsless_fields() -> None:
    density = np.asarray([0.5, 1.0, 1.5])
    blends = material_blend_arrays_for_density(density)
    discrete = discretize_material_blends(blends, mixture_bins=5)
    media_set = build_pegsless_media_set(discrete)

    rendered = render_media_definition(media_set)

    assert ":start media definition:" in rendered
    assert ":stop media definition:" in rendered
    assert "elements =" in rendered
    assert "mass fractions =" in rendered
    assert "bulk density =" in rendered
    assert "bremsstrahlung correction = NRC" in rendered
    assert "ivalue" not in rendered.lower()
