from pathlib import Path

import pytest
from pydicom.dataset import Dataset, FileDataset, FileMetaDataset
from pydicom.sequence import Sequence
from pydicom.uid import ExplicitVRLittleEndian, RTStructureSetStorage, generate_uid

from indep_mu.dicom.rtstruct import list_rtstruct_rois


def _write_rtstruct(path: Path) -> None:
    meta = FileMetaDataset()
    meta.MediaStorageSOPClassUID = RTStructureSetStorage
    meta.MediaStorageSOPInstanceUID = generate_uid()
    meta.TransferSyntaxUID = ExplicitVRLittleEndian

    ds = FileDataset(str(path), {}, file_meta=meta, preamble=b"\0" * 128)
    ds.SOPClassUID = RTStructureSetStorage
    ds.SOPInstanceUID = meta.MediaStorageSOPInstanceUID
    ds.Modality = "RTSTRUCT"

    frame_uid = generate_uid()
    body = Dataset()
    body.ROINumber = 1
    body.ROIName = "patient"
    body.ReferencedFrameOfReferenceUID = frame_uid

    ptv = Dataset()
    ptv.ROINumber = 2
    ptv.ROIName = "PTV_6996"
    ptv.ReferencedFrameOfReferenceUID = frame_uid

    ds.StructureSetROISequence = Sequence([body, ptv])
    ds.save_as(path, enforce_file_format=True)


def test_list_rtstruct_rois(tmp_path: Path) -> None:
    path = tmp_path / "struct.dcm"
    _write_rtstruct(path)

    rois = list_rtstruct_rois(path)

    assert [(item.roi_number, item.roi_name) for item in rois] == [
        (1, "patient"),
        (2, "PTV_6996"),
    ]


def test_duplicate_roi_numbers_are_rejected(tmp_path: Path) -> None:
    path = tmp_path / "struct.dcm"
    _write_rtstruct(path)

    import pydicom

    ds = pydicom.dcmread(path)
    ds.StructureSetROISequence[1].ROINumber = 1
    ds.save_as(path, enforce_file_format=True)

    with pytest.raises(ValueError, match="Duplicate ROINumber"):
        list_rtstruct_rois(path)
