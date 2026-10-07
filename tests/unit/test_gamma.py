import sys
from types import SimpleNamespace

import numpy as np
import pytest

from indep_mu.comparison3d.dose_difference import (
    DoseDifferenceResult,
    DoseDifferenceSummary,
)
from indep_mu.comparison3d.gamma import (
    GammaConfig,
    gamma_pass_rate,
    run_gamma_on_rtdose_grid,
)
from indep_mu.dicom.rtdose import DoseGeometry, RtDose


def _reference() -> RtDose:
    geometry = DoseGeometry(
        rows=2,
        columns=2,
        frames=2,
        pixel_spacing_mm=(2.0, 3.0),
        image_orientation_patient=(1.0, 0.0, 0.0, 0.0, 1.0, 0.0),
        image_position_patient_mm=(0.0, 0.0, 0.0),
        frame_offsets_mm=np.asarray([0.0, 2.5]),
        frame_positions_patient_mm=np.asarray(
            [[0.0, 0.0, 0.0], [0.0, 0.0, 2.5]]
        ),
        frame_offset_mode="RELATIVE",
    )
    dose = np.full((2, 2, 2), 10.0)
    return RtDose(
        dose=dose,
        geometry=geometry,
        dose_units="GY",
        dose_type="PHYSICAL",
        dose_summation_type="PLAN",
        frame_of_reference_uid="1.2.3",
        sop_instance_uid="1.2.3.4",
        referenced_rtplan_uids=(),
        referenced_beam_numbers=(),
    )


def _comparison(reference: RtDose) -> DoseDifferenceResult:
    evaluated = reference.dose_gy * 1.01
    difference = evaluated - reference.dose_gy
    relative = 100.0 * difference / reference.dose_gy
    mask = np.ones(reference.shape, dtype=bool)
    summary = DoseDifferenceSummary(
        evaluated_voxel_count=int(np.prod(reference.shape)),
        threshold_gy=1.0,
        mean_difference_gy=float(np.mean(difference)),
        mean_absolute_difference_gy=float(np.mean(np.abs(difference))),
        rms_difference_gy=float(np.sqrt(np.mean(difference**2))),
        max_absolute_difference_gy=float(np.max(np.abs(difference))),
        mean_relative_difference_percent=float(np.mean(relative)),
        p95_absolute_relative_difference_percent=float(
            np.percentile(np.abs(relative), 95)
        ),
    )
    return DoseDifferenceResult(
        reference_dose_gy=reference.dose_gy,
        evaluated_dose_gy=evaluated,
        difference_gy=difference,
        relative_difference_percent=relative,
        evaluation_mask=mask,
        summary=summary,
    )


def test_global_gamma_requires_explicit_normalisation() -> None:
    with pytest.raises(ValueError, match="explicit positive"):
        GammaConfig(
            dose_percent_threshold=3.0,
            distance_mm_threshold=3.0,
            lower_percent_dose_cutoff=10.0,
        )


def test_gamma_pass_rate_respects_mask() -> None:
    gamma = np.asarray(
        [
            [[0.5, 1.2], [0.8, np.nan]],
            [[1.0, 1.1], [0.2, 0.4]],
        ]
    )
    mask = np.zeros_like(gamma, dtype=bool)
    mask[0] = True

    analysed, passed, rate = gamma_pass_rate(gamma, mask=mask)

    assert analysed == 3
    assert passed == 2
    assert rate == pytest.approx(100.0 * 2.0 / 3.0)


def test_pymedphys_wrapper_receives_explicit_gamma_settings(monkeypatch) -> None:
    reference = _reference()
    comparison = _comparison(reference)
    calls = {}

    def fake_gamma(
        axes_reference,
        dose_reference,
        axes_evaluation,
        dose_evaluation,
        **kwargs,
    ):
        calls["axes_reference"] = axes_reference
        calls["dose_reference"] = dose_reference
        calls["axes_evaluation"] = axes_evaluation
        calls["dose_evaluation"] = dose_evaluation
        calls["kwargs"] = kwargs
        return np.full(reference.shape, 0.75)

    monkeypatch.setitem(
        sys.modules,
        "pymedphys",
        SimpleNamespace(gamma=fake_gamma),
    )

    config = GammaConfig(
        dose_percent_threshold=3.0,
        distance_mm_threshold=3.0,
        lower_percent_dose_cutoff=10.0,
        local_gamma=False,
        global_normalisation_gy=40.05,
        interp_fraction=10,
        max_gamma=2.0,
    )
    result = run_gamma_on_rtdose_grid(reference, comparison, config)

    assert result.pass_rate_percent == pytest.approx(100.0)
    assert calls["kwargs"]["dose_percent_threshold"] == pytest.approx(3.0)
    assert calls["kwargs"]["distance_mm_threshold"] == pytest.approx(3.0)
    assert calls["kwargs"]["lower_percent_dose_cutoff"] == pytest.approx(10.0)
    assert calls["kwargs"]["global_normalisation"] == pytest.approx(40.05)
    assert calls["kwargs"]["local_gamma"] is False
    assert calls["kwargs"]["skip_once_passed"] is True

    z, y, x = calls["axes_reference"]
    np.testing.assert_allclose(z, [0.0, 2.5])
    np.testing.assert_allclose(y, [0.0, 2.0])
    np.testing.assert_allclose(x, [0.0, 3.0])
