from __future__ import annotations

from pathlib import Path

from pydicom.dataset import Dataset, FileDataset, FileMetaDataset
from pydicom.sequence import Sequence
from pydicom.uid import (
    CTImageStorage,
    ExplicitVRLittleEndian,
    RTDoseStorage,
    RTPlanStorage,
    RTStructureSetStorage,
    generate_uid,
)

from indep_mu.dicom.case import discover_dicom_case


def _file_dataset(path: Path, sop_class_uid: str, modality: str) -> FileDataset:
    meta = FileMetaDataset()
    meta.MediaStorageSOPClassUID = sop_class_uid
    meta.MediaStorageSOPInstanceUID = generate_uid()
    meta.TransferSyntaxUID = ExplicitVRLittleEndian

    ds = FileDataset(str(path), {}, file_meta=meta, preamble=b"\0" * 128)
    ds.SOPClassUID = sop_class_uid
    ds.SOPInstanceUID = meta.MediaStorageSOPInstanceUID
    ds.Modality = modality
    return ds


def _write_ct(
    path: Path,
    *,
    series_uid: str,
    frame_uid: str,
) -> str:
    ds = _file_dataset(path, CTImageStorage, "CT")
    ds.StudyInstanceUID = generate_uid()
    ds.SeriesInstanceUID = series_uid
    ds.FrameOfReferenceUID = frame_uid
    ds.save_as(path, enforce_file_format=True)
    return str(ds.SOPInstanceUID)


def _write_struct(
    path: Path,
    *,
    frame_uid: str,
    ct_series_uid: str,
    ct_sop_uid: str,
) -> str:
    ds = _file_dataset(path, RTStructureSetStorage, "RTSTRUCT")
    ds.StudyInstanceUID = generate_uid()
    ds.SeriesInstanceUID = generate_uid()

    frame = Dataset()
    frame.FrameOfReferenceUID = frame_uid
    study = Dataset()
    series = Dataset()
    series.SeriesInstanceUID = ct_series_uid
    image = Dataset()
    image.ReferencedSOPClassUID = CTImageStorage
    image.ReferencedSOPInstanceUID = ct_sop_uid
    series.ContourImageSequence = Sequence([image])
    study.RTReferencedSeriesSequence = Sequence([series])
    frame.RTReferencedStudySequence = Sequence([study])
    ds.ReferencedFrameOfReferenceSequence = Sequence([frame])

    roi = Dataset()
    roi.ROINumber = 1
    roi.ROIName = "patient"
    roi.ReferencedFrameOfReferenceUID = frame_uid
    ds.StructureSetROISequence = Sequence([roi])

    contour_image = Dataset()
    contour_image.ReferencedSOPClassUID = CTImageStorage
    contour_image.ReferencedSOPInstanceUID = ct_sop_uid
    contour = Dataset()
    contour.ContourImageSequence = Sequence([contour_image])
    roi_contour = Dataset()
    roi_contour.ReferencedROINumber = 1
    roi_contour.ContourSequence = Sequence([contour])
    ds.ROIContourSequence = Sequence([roi_contour])

    ds.save_as(path, enforce_file_format=True)
    return str(ds.SOPInstanceUID)


def _write_plan(path: Path, *, struct_uid: str) -> str:
    ds = _file_dataset(path, RTPlanStorage, "RTPLAN")
    ds.StudyInstanceUID = generate_uid()
    ds.SeriesInstanceUID = generate_uid()

    ref = Dataset()
    ref.ReferencedSOPClassUID = RTStructureSetStorage
    ref.ReferencedSOPInstanceUID = struct_uid
    ds.ReferencedStructureSetSequence = Sequence([ref])

    ds.save_as(path, enforce_file_format=True)
    return str(ds.SOPInstanceUID)


def _write_dose(
    path: Path,
    *,
    frame_uid: str,
    plan_uid: str,
) -> None:
    ds = _file_dataset(path, RTDoseStorage, "RTDOSE")
    ds.StudyInstanceUID = generate_uid()
    ds.SeriesInstanceUID = generate_uid()
    ds.FrameOfReferenceUID = frame_uid

    ref = Dataset()
    ref.ReferencedSOPClassUID = RTPlanStorage
    ref.ReferencedSOPInstanceUID = plan_uid
    ds.ReferencedRTPlanSequence = Sequence([ref])

    ds.save_as(path, enforce_file_format=True)


def test_discovers_referenced_ct_not_unrelated_series(tmp_path: Path) -> None:
    frame_uid = generate_uid()
    selected_series = generate_uid()
    unrelated_series = generate_uid()

    selected_sop = _write_ct(
        tmp_path / "ct_selected.dcm",
        series_uid=selected_series,
        frame_uid=frame_uid,
    )
    _write_ct(
        tmp_path / "ct_unrelated.dcm",
        series_uid=unrelated_series,
        frame_uid=generate_uid(),
    )

    struct_uid = _write_struct(
        tmp_path / "struct.dcm",
        frame_uid=frame_uid,
        ct_series_uid=selected_series,
        ct_sop_uid=selected_sop,
    )
    plan_uid = _write_plan(tmp_path / "plan.dcm", struct_uid=struct_uid)
    _write_dose(
        tmp_path / "dose.dcm",
        frame_uid=frame_uid,
        plan_uid=plan_uid,
    )

    case = discover_dicom_case(tmp_path)

    assert case.ct_series.series_instance_uid == selected_series
    assert len(case.ct_series.instances) == 1
    assert len(case.rtdoses) == 1
    assert case.passed


def test_flags_rtdose_frame_mismatch(tmp_path: Path) -> None:
    frame_uid = generate_uid()
    selected_series = generate_uid()

    selected_sop = _write_ct(
        tmp_path / "ct.dcm",
        series_uid=selected_series,
        frame_uid=frame_uid,
    )
    struct_uid = _write_struct(
        tmp_path / "struct.dcm",
        frame_uid=frame_uid,
        ct_series_uid=selected_series,
        ct_sop_uid=selected_sop,
    )
    plan_uid = _write_plan(tmp_path / "plan.dcm", struct_uid=struct_uid)
    _write_dose(
        tmp_path / "dose.dcm",
        frame_uid=generate_uid(),
        plan_uid=plan_uid,
    )

    case = discover_dicom_case(tmp_path)

    assert not case.passed
    assert any(item.code == "RTDOSE_FRAME_MISMATCH" for item in case.errors)
