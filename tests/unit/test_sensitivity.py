import numpy as np

from indep_mu.patient_model.materials import material_blend_arrays_for_density
from indep_mu.patient_model.sensitivity import run_discretization_sensitivity


def test_more_mixture_bins_reduce_or_preserve_fraction_error() -> None:
    density = np.linspace(0.01, 2.5, 500)
    blends = material_blend_arrays_for_density(density)

    results = run_discretization_sensitivity(
        blends,
        density,
        candidate_mixture_bins=(5, 9, 19),
    )

    assert all(result.feasible for result in results)
    assert results[1].max_abs_fraction_error <= results[0].max_abs_fraction_error
    assert results[2].max_abs_fraction_error <= results[1].max_abs_fraction_error


def test_sensitivity_reports_density_scaling() -> None:
    density = np.asarray([0.001, 0.2, 0.9, 1.0, 1.5, 2.5])
    blends = material_blend_arrays_for_density(density)

    result = run_discretization_sensitivity(
        blends,
        density,
        candidate_mixture_bins=(5,),
    )[0]

    assert result.feasible
    assert result.media_count is not None
    assert result.rhof_min is not None
    assert result.rhof_max is not None
    assert result.rhof_min > 0
    assert result.rhof_max > 0
