import pytest

from indep_mu.beam_model.tps_reference import (
    TpsReferenceDataset,
    TpsReferenceDose,
    TpsReferenceGeometry,
)
from indep_mu.beam_model.reference_check import (
    compare_independent_dose_to_tps_reference,
)
from indep_mu.beam_model.commissioning import FieldSize


def _reference() -> TpsReferenceDataset:
    return TpsReferenceDataset(
        machine_model_id="VERSA_HD",
        tps_name="Elekta Monaco",
        tps_version="6.1.4",
        geometry=TpsReferenceGeometry(
            field=FieldSize(100.0, 100.0),
            ssd_mm=900.0,
            depth_mm=100.0,
            delivered_mu=100.0,
        ),
        energies=(TpsReferenceDose("6MV", 0.995),),
    )


def test_reference_check_does_not_renormalize() -> None:
    result = compare_independent_dose_to_tps_reference(
        energy_model_id="6MV",
        independent_dose_gy=1.000,
        tps_reference=_reference(),
    )

    assert result.independent_dose_gy == pytest.approx(1.000)
    assert result.tps_reference_dose_gy == pytest.approx(0.995)
    assert result.difference_percent == pytest.approx((1.000 / 0.995 - 1.0) * 100.0)


def test_reference_check_rejects_nonpositive_independent_dose() -> None:
    with pytest.raises(ValueError):
        compare_independent_dose_to_tps_reference(
            energy_model_id="6MV",
            independent_dose_gy=0.0,
            tps_reference=_reference(),
        )
