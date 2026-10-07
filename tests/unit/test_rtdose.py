from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from pydicom.dataset import FileDataset, FileMetaDataset
from pydicom.uid import ExplicitVRLittleEndian, RTDoseStorage, generate_uid

from indep_mu.dicom.rtdose import load_rtdose


def _write_dose(
    path: Path,
    *,
    offsets: list[float],
    image_position=(10.0, 20.0, 30.0),
    orientation=(1.0, 0.0, 0.0, 0.0, 1.0, 0.0),
) -> None:
    meta = FileMetaDataset()
    meta.MediaStorageSOPClassUID = RTDoseStorage
    meta.MediaStorageSOPInstanceUID = generate_uid()
    meta.TransferSyntaxUID = ExplicitVRLittleEndian

    ds = FileDataset(str(path), {}, file_meta=meta, preamble=b"\0" * 128)
    ds.SOPClassUID = RTDoseStorage
    ds.SOPInstanceUID = meta.MediaStorageSOPInstanceUID
    ds.Modality = "RTDOSE"
    ds.FrameOfReferenceUID = generate_uid()

    ds.Rows = 2
    ds.Columns = 3
    ds.NumberOfFrames = len(offsets)
    ds.PixelSpacing = [2.0, 3.0]
    ds.ImageOrientationPatient = list(orientation)
    ds.ImagePositionPatient = list(image_position)
    ds.GridFrameOffsetVector = offsets

    ds.DoseUnits = "GY"
    ds.DoseType = "PHYSICAL"
    ds.DoseSummationType = "PLAN"
    ds.DoseGridScaling = 0.01

    ds.SamplesPerPixel = 1
    ds.PhotometricInterpretation = "MONOCHROME2"
    ds.BitsAllocated = 16
    ds.BitsStored = 16
    ds.HighBit = 15
    ds.PixelRepresentation = 0

    pixels = np.arange(
        len(offsets) * ds.Rows * ds.Columns,
        dtype=np.uint16,
    ).reshape((len(offsets), ds.Rows, ds.Columns))
    ds.PixelData = pixels.tobytes()
    ds.save_as(path, enforce_file_format=True)


def test_relative_grid_frame_offset_vector(tmp_path: Path) -> None:
    path = tmp_path / "dose.dcm"
    _write_dose(path, offsets=[0.0, 2.5])

    dose = load_rtdose(path)

    assert dose.shape == (2, 2, 3)
    assert dose.geometry.frame_offset_mode == "RELATIVE"
    np.testing.assert_allclose(
        dose.geometry.frame_positions_patient_mm,
        [[10.0, 20.0, 30.0], [10.0, 20.0, 32.5]],
    )
    np.testing.assert_allclose(
        dose.geometry.voxel_center_patient_mm(1, 1, 2),
        [16.0, 22.0, 32.5],
    )
    assert dose.dose[1, 1, 2] == pytest.approx(0.11)


def test_legacy_absolute_patient_z_offsets(tmp_path: Path) -> None:
    path = tmp_path / "dose.dcm"
    _write_dose(path, offsets=[30.0, 32.5])

    dose = load_rtdose(path)

    assert dose.geometry.frame_offset_mode == "ABSOLUTE_PATIENT_Z_LEGACY"
    np.testing.assert_allclose(
        dose.geometry.frame_positions_patient_mm[:, 2],
        [30.0, 32.5],
    )


def test_nonzero_relative_offset_is_rejected_for_oblique_grid(
    tmp_path: Path,
) -> None:
    path = tmp_path / "dose.dcm"
    _write_dose(
        path,
        offsets=[30.0, 32.5],
        orientation=(0.0, 1.0, 0.0, 0.0, 0.0, 1.0),
    )

    with pytest.raises(ValueError, match="neither DICOM relative option"):
        load_rtdose(path)
