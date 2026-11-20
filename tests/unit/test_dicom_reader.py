"""Regression coverage for Monaco headerless implicit-VR DICOM exports."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pydicom
import pytest
from pydicom.dataset import Dataset, FileDataset, FileMetaDataset
from pydicom.uid import (
    CTImageStorage,
    ExplicitVRLittleEndian,
    generate_uid,
)

from indep_mu.dicom.reader import dcmread_rt


def _ct_without_part10(path: Path) -> None:
    ds = Dataset()
    ds.SOPClassUID = CTImageStorage
    ds.SOPInstanceUID = generate_uid()
    ds.Modality = "CT"
    ds.Rows = 2
    ds.Columns = 3
    ds.SamplesPerPixel = 1
    ds.PhotometricInterpretation = "MONOCHROME2"
    ds.BitsAllocated = 16
    ds.BitsStored = 16
    ds.HighBit = 15
    ds.PixelRepresentation = 1
    ds.RescaleSlope = 1
    ds.RescaleIntercept = -1024
    ds.PixelData = np.arange(6, dtype="<i2").tobytes()

    # No preamble, no DICM marker, no File Meta Information, implicit VR.
    pydicom.dcmwrite(
        path,
        ds,
        implicit_vr=True,
        little_endian=True,
        enforce_file_format=False,
    )


def test_headerless_monaco_ct_is_readable_with_pixel_data(tmp_path: Path) -> None:
    path = tmp_path / "raw_legacy.dcm"
    _ct_without_part10(path)
    data = path.read_bytes()
    assert data[128:132] != b"DICM"
    assert data[:2] == b"\x08\x00"

    header = dcmread_rt(path, stop_before_pixels=True)
    assert header.Modality == "CT"
    assert not hasattr(header, "PixelData")

    dataset = dcmread_rt(path)
    assert dataset.file_meta.TransferSyntaxUID.is_implicit_VR
    np.testing.assert_array_equal(
        dataset.pixel_array,
        np.arange(6, dtype=np.int16).reshape(2, 3),
    )


def test_part10_ct_uses_normal_strict_path(tmp_path: Path) -> None:
    path = tmp_path / "standard.dcm"
    meta = FileMetaDataset()
    meta.MediaStorageSOPClassUID = CTImageStorage
    meta.MediaStorageSOPInstanceUID = generate_uid()
    meta.TransferSyntaxUID = ExplicitVRLittleEndian

    ds = FileDataset(str(path), {}, file_meta=meta, preamble=b"\0" * 128)
    ds.SOPClassUID = CTImageStorage
    ds.SOPInstanceUID = meta.MediaStorageSOPInstanceUID
    ds.Modality = "CT"
    ds.save_as(path, enforce_file_format=True)

    restored = dcmread_rt(path)
    assert restored.SOPInstanceUID == ds.SOPInstanceUID
    assert restored.file_meta.TransferSyntaxUID == ExplicitVRLittleEndian


def test_random_binary_is_not_accepted_by_force_fallback(tmp_path: Path) -> None:
    path = tmp_path / "garbage.bin"
    path.write_bytes(b"NOT A DICOM FILE" * 64)
    with pytest.raises(ValueError, match="does not resemble"):
        dcmread_rt(path)


def test_raw_dataset_with_unknown_class_is_rejected(tmp_path: Path) -> None:
    path = tmp_path / "unsupported.dcm"
    ds = Dataset()
    ds.SOPClassUID = generate_uid()
    ds.SOPInstanceUID = generate_uid()
    ds.Modality = "CT"
    pydicom.dcmwrite(
        path,
        ds,
        implicit_vr=True,
        little_endian=True,
        enforce_file_format=False,
    )
    with pytest.raises(ValueError, match="approved CT/RT"):
        dcmread_rt(path)
