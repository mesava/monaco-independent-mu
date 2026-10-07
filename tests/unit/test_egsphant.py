from pathlib import Path

import numpy as np

from indep_mu.dicom.ct import CtGeometry, CtSeries
from indep_mu.patient_model.discretize import discretize_material_blends
from indep_mu.patient_model.egsphant import (
    EGSPHANT_ENCODING,
    build_egsphant_geometry,
    read_egsphant,
    write_egsphant,
)
from indep_mu.patient_model.materials import material_blend_arrays_for_density
from indep_mu.patient_model.model import CropBounds, PatientModel


def _synthetic_ct_and_model() -> tuple[CtSeries, PatientModel]:
    geometry = CtGeometry(
        rows=2,
        columns=3,
        pixel_spacing_mm=(2.0, 1.0),
        slice_spacing_mm=2.5,
        image_orientation_patient=(1.0, 0.0, 0.0, 0.0, 1.0, 0.0),
        first_image_position_patient_mm=(10.0, 20.0, 30.0),
        last_image_position_patient_mm=(10.0, 20.0, 32.5),
        patient_position="HFS",
        frame_of_reference_uid="1.2.3",
        series_instance_uid="1.2.3.4",
        kvp=120.0,
    )
    hu = np.zeros((2, 2, 3), dtype=np.float32)
    ct = CtSeries(
        hu=hu,
        geometry=geometry,
        slice_positions_mm=np.asarray([30.0, 32.5]),
        image_positions_patient_mm=np.asarray(
            [[10.0, 20.0, 30.0], [10.0, 20.0, 32.5]]
        ),
        sop_instance_uids=("1", "2"),
    )

    density = np.asarray(
        [
            [[0.001, 0.50, 0.92], [1.00, 1.20, 1.85]],
            [[0.001, 0.70, 0.95], [1.04, 1.50, 2.57]],
        ],
        dtype=np.float32,
    )
    blends = material_blend_arrays_for_density(density)
    model = PatientModel(
        hu=hu,
        red=np.ones_like(hu),
        mass_density_g_cm3=density,
        patient_mask=np.ones_like(hu, dtype=bool),
        out_of_range_code=np.zeros_like(hu, dtype=np.int8),
        material_blends=blends,
        crop_bounds=CropBounds(0, 2, 0, 2, 0, 3),
        voxel_volume_mm3=5.0,
    )
    return ct, model


def test_encoding_matches_egsnrc_macro() -> None:
    assert EGSPHANT_ENCODING == "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"
    assert EGSPHANT_ENCODING[1] == "1"
    assert EGSPHANT_ENCODING[61] == "z"


def test_local_geometry_preserves_dicom_spacing() -> None:
    ct, model = _synthetic_ct_and_model()
    geo = build_egsphant_geometry(ct, model)

    np.testing.assert_allclose(np.diff(geo.x_bound_cm), 0.1)
    np.testing.assert_allclose(np.diff(geo.y_bound_cm), 0.2)
    np.testing.assert_allclose(np.diff(geo.z_bound_cm), 0.25)
    np.testing.assert_allclose(geo.x_bound_cm[[0, -1]], [0.95, 1.25])
    np.testing.assert_allclose(geo.y_bound_cm[[0, -1]], [1.9, 2.3])
    np.testing.assert_allclose(geo.z_bound_cm[[0, -1]], [2.875, 3.375])


def test_egsphant_round_trip(tmp_path: Path) -> None:
    ct, model = _synthetic_ct_and_model()
    discrete = discretize_material_blends(model.material_blends, mixture_bins=5)
    path = tmp_path / "test.egsphant"

    written_geometry = write_egsphant(
        path,
        ct=ct,
        model=model,
        discrete=discrete,
    )
    parsed = read_egsphant(path)

    assert parsed.media_names == tuple(item.egsphant_name for item in discrete.media)
    np.testing.assert_array_equal(parsed.medium_index, discrete.medium_index)
    np.testing.assert_allclose(
        parsed.mass_density_g_cm3,
        model.mass_density_g_cm3,
        rtol=0.0,
        atol=1e-7,
    )
    np.testing.assert_allclose(parsed.geometry.x_bound_cm, written_geometry.x_bound_cm)
    np.testing.assert_allclose(parsed.geometry.y_bound_cm, written_geometry.y_bound_cm)
    np.testing.assert_allclose(parsed.geometry.z_bound_cm, written_geometry.z_bound_cm)
