import numpy as np
import pytest

from indep_mu.dicom.ct import CtGeometry, CtSeries
from indep_mu.patient_model.external_mask import (
    derive_external_mask_from_ct,
    patient_point_to_ct_index_zyx,
)


def _ct() -> CtSeries:
    shape = (7, 11, 11)
    hu = np.full(shape, -1000, dtype=np.int16)

    zz, yy, xx = np.indices(shape)
    body = (
        ((zz - 3.0) / 2.5) ** 2
        + ((yy - 5.0) / 4.0) ** 2
        + ((xx - 5.0) / 4.0) ** 2
        <= 1.0
    )
    hu[body] = 0

    # Enclosed low-density cavity: should still be inside the external mask.
    hu[3, 4:7, 4:7] = -1000

    # Detached synthetic couch/table object.
    hu[1:6, 10, 1:10] = 0

    geometry = CtGeometry(
        rows=11,
        columns=11,
        pixel_spacing_mm=(1.0, 1.0),
        slice_spacing_mm=1.0,
        image_orientation_patient=(1.0, 0.0, 0.0, 0.0, 1.0, 0.0),
        first_image_position_patient_mm=(0.0, 0.0, 0.0),
        last_image_position_patient_mm=(0.0, 0.0, 6.0),
        patient_position="HFS",
        frame_of_reference_uid="1.2.3",
        series_instance_uid="1.2.3.4",
        kvp=120.0,
    )
    positions = np.asarray(
        [[0.0, 0.0, float(z)] for z in range(shape[0])],
        dtype=np.float64,
    )
    return CtSeries(
        hu=hu,
        geometry=geometry,
        slice_positions_mm=np.arange(shape[0], dtype=np.float64),
        image_positions_patient_mm=positions,
        sop_instance_uids=tuple(str(index) for index in range(shape[0])),
    )


def test_patient_point_maps_to_ct_voxel() -> None:
    ct = _ct()
    assert patient_point_to_ct_index_zyx(ct, (5.0, 5.0, 3.0)) == (3, 5, 5)


def test_seeded_external_mask_excludes_detached_couch_and_fills_hole() -> None:
    ct = _ct()

    result = derive_external_mask_from_ct(
        ct,
        seed_patient_mm=(5.0, 5.0, 3.0),
        threshold_hu=-500.0,
        closing_iterations=0,
        fill_holes_per_slice=True,
    )

    # Central cavity becomes part of the external patient envelope.
    assert result.mask[3, 5, 5]

    # Detached couch remains excluded.
    assert not np.any(result.mask[:, 10, :])

    assert result.diagnostics.connected_component_count >= 2
    assert not result.diagnostics.touches_ct_border


def test_seed_in_air_is_rejected() -> None:
    ct = _ct()

    with pytest.raises(ValueError, match="not inside"):
        derive_external_mask_from_ct(
            ct,
            seed_patient_mm=(0.0, 0.0, 0.0),
            threshold_hu=-500.0,
            closing_iterations=0,
        )
