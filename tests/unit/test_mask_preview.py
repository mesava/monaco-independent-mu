from pathlib import Path

import numpy as np
import pytest

from indep_mu.dicom.ct import CtGeometry, CtSeries
from indep_mu.patient_model.mask_preview import render_external_mask_montage


def _synthetic_ct() -> tuple[CtSeries, np.ndarray]:
    hu = np.full((5, 12, 12), -1000.0, dtype=np.float32)
    mask = np.zeros_like(hu, dtype=bool)
    mask[1:4, 3:9, 2:10] = True
    hu[mask] = 0.0

    geometry = CtGeometry(
        rows=12,
        columns=12,
        pixel_spacing_mm=(2.0, 2.0),
        slice_spacing_mm=2.5,
        image_orientation_patient=(1.0, 0.0, 0.0, 0.0, 1.0, 0.0),
        first_image_position_patient_mm=(0.0, 0.0, 0.0),
        last_image_position_patient_mm=(0.0, 0.0, 10.0),
        patient_position="HFS",
        frame_of_reference_uid="1.2.3",
        series_instance_uid="1.2.3.4",
        kvp=120.0,
    )
    ct = CtSeries(
        hu=hu,
        geometry=geometry,
        slice_positions_mm=np.arange(5, dtype=float) * 2.5,
        image_positions_patient_mm=np.column_stack(
            (np.zeros(5), np.zeros(5), np.arange(5, dtype=float) * 2.5)
        ),
        sop_instance_uids=tuple(str(i) for i in range(5)),
    )
    return ct, mask


def test_local_qc_png_does_not_include_dicom_identifiers(tmp_path: Path) -> None:
    pytest.importorskip("matplotlib")
    ct, mask = _synthetic_ct()
    dest = tmp_path / "local_mask_preview.png"

    indices = render_external_mask_montage(ct, mask, dest)

    assert indices == (1, 2, 3)
    assert dest.read_bytes().startswith(b"\x89PNG\r\n\x1a\n")
    assert len(dest.read_bytes()) > 1_000


def test_qc_png_does_not_overwrite_silently(tmp_path: Path) -> None:
    ct, mask = _synthetic_ct()
    dest = tmp_path / "existing.png"
    dest.write_bytes(b"keep")

    with pytest.raises(FileExistsError):
        render_external_mask_montage(ct, mask, dest)

    assert dest.read_bytes() == b"keep"


def test_qc_png_rejects_empty_mask(tmp_path: Path) -> None:
    ct, mask = _synthetic_ct()
    with pytest.raises(ValueError, match="empty patient mask"):
        render_external_mask_montage(
            ct,
            np.zeros_like(mask),
            tmp_path / "empty.png",
        )



def test_local_png_prioritizes_suspicious_axial_slices(tmp_path: Path) -> None:
    pytest.importorskip("matplotlib")
    ct, mask = _synthetic_ct()
    dest = tmp_path / "priority_mask_preview.png"

    indices = render_external_mask_montage(
        ct,
        mask,
        dest,
        panel_count=2,
        priority_slice_indices=(3,),
    )

    assert indices == (1, 3)
    assert dest.stat().st_size > 1000
