from __future__ import annotations

from pathlib import Path

import numpy as np
import pydicom
import pytest
from pydicom.dataset import FileDataset, FileMetaDataset
from pydicom.uid import CTImageStorage, ExplicitVRLittleEndian, generate_uid

from indep_mu.dicom.ct import load_ct_series


def _write_ct_slice(
    path: Path,
    *,
    series_uid: str,
    frame_uid: str,
    z_mm: float,
    value: int,
    slope: float = 1.0,
    intercept: float = -1024.0,
) -> None:
    meta = FileMetaDataset()
    meta.MediaStorageSOPClassUID = CTImageStorage
    meta.MediaStorageSOPInstanceUID = generate_uid()
    meta.TransferSyntaxUID = ExplicitVRLittleEndian

    ds = FileDataset(str(path), {}, file_meta=meta, preamble=b"\0" * 128)
    ds.SOPClassUID = CTImageStorage
    ds.SOPInstanceUID = meta.MediaStorageSOPInstanceUID
    ds.Modality = "CT"
    ds.SeriesInstanceUID = series_uid
    ds.StudyInstanceUID = generate_uid()
    ds.FrameOfReferenceUID = frame_uid

    ds.Rows = 4
    ds.Columns = 5
    ds.PixelSpacing = [1.0, 1.0]
    ds.ImageOrientationPatient = [1, 0, 0, 0, 1, 0]
    ds.ImagePositionPatient = [0, 0, z_mm]
    ds.PatientPosition = "HFS"
    ds.KVP = 120

    ds.SamplesPerPixel = 1
    ds.PhotometricInterpretation = "MONOCHROME2"
    ds.BitsAllocated = 16
    ds.BitsStored = 16
    ds.HighBit = 15
    ds.PixelRepresentation = 1
    ds.RescaleSlope = slope
    ds.RescaleIntercept = intercept

    pixels = np.full((ds.Rows, ds.Columns), value, dtype=np.int16)
    ds.PixelData = pixels.tobytes()
    ds.save_as(path, enforce_file_format=True)


def test_ct_series_is_sorted_and_converted_to_hu(tmp_path: Path) -> None:
    series_uid = generate_uid()
    frame_uid = generate_uid()

    _write_ct_slice(
        tmp_path / "slice2.dcm",
        series_uid=series_uid,
        frame_uid=frame_uid,
        z_mm=2.5,
        value=1100,
    )
    _write_ct_slice(
        tmp_path / "slice0.dcm",
        series_uid=series_uid,
        frame_uid=frame_uid,
        z_mm=-2.5,
        value=1000,
    )
    _write_ct_slice(
        tmp_path / "slice1.dcm",
        series_uid=series_uid,
        frame_uid=frame_uid,
        z_mm=0.0,
        value=1024,
    )

    ct = load_ct_series(tmp_path)

    assert ct.hu.shape == (3, 4, 5)
    np.testing.assert_allclose(ct.slice_positions_mm, [-2.5, 0.0, 2.5])
    assert ct.geometry.slice_spacing_mm == pytest.approx(2.5)
    assert ct.geometry.kvp == pytest.approx(120.0)
    assert float(ct.hu[0, 0, 0]) == pytest.approx(-24.0)
    assert float(ct.hu[1, 0, 0]) == pytest.approx(0.0)
    assert float(ct.hu[2, 0, 0]) == pytest.approx(76.0)


def test_mixed_frame_of_reference_is_rejected(tmp_path: Path) -> None:
    series_uid = generate_uid()

    _write_ct_slice(
        tmp_path / "a.dcm",
        series_uid=series_uid,
        frame_uid=generate_uid(),
        z_mm=0.0,
        value=1024,
    )
    _write_ct_slice(
        tmp_path / "b.dcm",
        series_uid=series_uid,
        frame_uid=generate_uid(),
        z_mm=2.5,
        value=1024,
    )

    with pytest.raises(ValueError, match="FrameOfReferenceUID"):
        load_ct_series(tmp_path)


def test_nonuniform_slice_spacing_is_rejected(tmp_path: Path) -> None:
    series_uid = generate_uid()
    frame_uid = generate_uid()

    for index, z in enumerate((0.0, 2.5, 5.2)):
        _write_ct_slice(
            tmp_path / f"{index}.dcm",
            series_uid=series_uid,
            frame_uid=frame_uid,
            z_mm=z,
            value=1024,
        )

    with pytest.raises(ValueError, match="Non-uniform CT slice spacing"):
        load_ct_series(tmp_path)



def test_requested_series_ignores_unrelated_ct_series(tmp_path: Path) -> None:
    selected_series = generate_uid()
    selected_frame = generate_uid()
    unrelated_series = generate_uid()

    for index, z in enumerate((0.0, 2.5, 5.0)):
        _write_ct_slice(
            tmp_path / f"selected_{index}.dcm",
            series_uid=selected_series,
            frame_uid=selected_frame,
            z_mm=z,
            value=1024 + index,
        )

    for index, z in enumerate((0.0, 5.0)):
        _write_ct_slice(
            tmp_path / f"unrelated_{index}.dcm",
            series_uid=unrelated_series,
            frame_uid=generate_uid(),
            z_mm=z,
            value=900,
        )

    ct = load_ct_series(
        tmp_path,
        series_instance_uid=selected_series,
    )

    assert ct.geometry.series_instance_uid == selected_series
    assert ct.geometry.frame_of_reference_uid == selected_frame
    assert ct.hu.shape[0] == 3


def test_mixed_series_without_explicit_selection_still_fails(tmp_path: Path) -> None:
    for series_index in range(2):
        series_uid = generate_uid()
        frame_uid = generate_uid()
        for slice_index, z in enumerate((0.0, 2.5)):
            _write_ct_slice(
                tmp_path / f"s{series_index}_{slice_index}.dcm",
                series_uid=series_uid,
                frame_uid=frame_uid,
                z_mm=z,
                value=1024,
            )

    with pytest.raises(ValueError, match="SeriesInstanceUID"):
        load_ct_series(tmp_path)
