from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from pydicom.dataset import Dataset, FileDataset, FileMetaDataset
from pydicom.sequence import Sequence
from pydicom.uid import (
    CTImageStorage,
    ExplicitVRLittleEndian,
    RTPlanStorage,
    RTStructureSetStorage,
    generate_uid,
)

from indep_mu.patient_model.egsphant import read_egsphant
from indep_mu.workflow.patient_build import build_patient_artifacts


def _base_dataset(path: Path, sop_class_uid: str, modality: str) -> FileDataset:
    meta = FileMetaDataset()
    meta.MediaStorageSOPClassUID = sop_class_uid
    meta.MediaStorageSOPInstanceUID = generate_uid()
    meta.TransferSyntaxUID = ExplicitVRLittleEndian

    ds = FileDataset(str(path), {}, file_meta=meta, preamble=b"\0" * 128)
    ds.SOPClassUID = sop_class_uid
    ds.SOPInstanceUID = meta.MediaStorageSOPInstanceUID
    ds.Modality = modality
    return ds


def _write_ct_slice(
    path: Path,
    *,
    series_uid: str,
    study_uid: str,
    frame_uid: str,
    z_mm: float,
) -> str:
    ds = _base_dataset(path, CTImageStorage, "CT")
    ds.SeriesInstanceUID = series_uid
    ds.StudyInstanceUID = study_uid
    ds.FrameOfReferenceUID = frame_uid
    ds.Rows = 4
    ds.Columns = 4
    ds.PixelSpacing = [1.0, 1.0]
    ds.ImageOrientationPatient = [1, 0, 0, 0, 1, 0]
    ds.ImagePositionPatient = [0.0, 0.0, z_mm]
    ds.PatientPosition = "HFS"
    ds.KVP = 120
    ds.RescaleSlope = 1.0
    ds.RescaleIntercept = -1024.0
    ds.SamplesPerPixel = 1
    ds.PhotometricInterpretation = "MONOCHROME2"
    ds.BitsAllocated = 16
    ds.BitsStored = 16
    ds.HighBit = 15
    ds.PixelRepresentation = 1
    ds.PixelData = np.full((4, 4), 1024, dtype=np.int16).tobytes()
    ds.save_as(path, enforce_file_format=True)
    return str(ds.SOPInstanceUID)


def _write_rtstruct(
    path: Path,
    *,
    study_uid: str,
    frame_uid: str,
    ct_series_uid: str,
    ct_sops: tuple[str, ...],
    z_values: tuple[float, ...],
) -> str:
    ds = _base_dataset(path, RTStructureSetStorage, "RTSTRUCT")
    ds.StudyInstanceUID = study_uid
    ds.SeriesInstanceUID = generate_uid()

    frame = Dataset()
    frame.FrameOfReferenceUID = frame_uid
    study = Dataset()
    study.ReferencedSOPInstanceUID = study_uid
    series = Dataset()
    series.SeriesInstanceUID = ct_series_uid
    series.ContourImageSequence = Sequence([])
    study.RTReferencedSeriesSequence = Sequence([series])
    frame.RTReferencedStudySequence = Sequence([study])
    ds.ReferencedFrameOfReferenceSequence = Sequence([frame])

    roi = Dataset()
    roi.ROINumber = 1
    roi.ROIName = "patient"
    roi.ReferencedFrameOfReferenceUID = frame_uid
    ds.StructureSetROISequence = Sequence([roi])

    contours = []
    for sop_uid, z_mm in zip(ct_sops, z_values):
        image = Dataset()
        image.ReferencedSOPClassUID = CTImageStorage
        image.ReferencedSOPInstanceUID = sop_uid

        contour = Dataset()
        contour.ContourGeometricType = "CLOSED_PLANAR"
        contour.NumberOfContourPoints = 4
        contour.ContourImageSequence = Sequence([image])
        contour.ContourData = [
            -0.5,
            -0.5,
            z_mm,
            3.5,
            -0.5,
            z_mm,
            3.5,
            3.5,
            z_mm,
            -0.5,
            3.5,
            z_mm,
        ]
        contours.append(contour)

    roi_contour = Dataset()
    roi_contour.ReferencedROINumber = 1
    roi_contour.ContourSequence = Sequence(contours)
    ds.ROIContourSequence = Sequence([roi_contour])

    ds.save_as(path, enforce_file_format=True)
    return str(ds.SOPInstanceUID)


def _write_rtplan(path: Path, *, study_uid: str, struct_uid: str) -> None:
    ds = _base_dataset(path, RTPlanStorage, "RTPLAN")
    ds.StudyInstanceUID = study_uid
    ds.SeriesInstanceUID = generate_uid()
    ds.RTPlanLabel = "TEST"
    ds.RTPlanName = "TEST"

    ref = Dataset()
    ref.ReferencedSOPClassUID = RTStructureSetStorage
    ref.ReferencedSOPInstanceUID = struct_uid
    ds.ReferencedStructureSetSequence = Sequence([ref])

    device = Dataset()
    device.RTBeamLimitingDeviceType = "ASYMY"
    device.NumberOfLeafJawPairs = 1
    device.SourceToBeamLimitingDeviceDistance = 470.0

    cp0_position = Dataset()
    cp0_position.RTBeamLimitingDeviceType = "ASYMY"
    cp0_position.LeafJawPositions = [-20.0, 20.0]

    cp0 = Dataset()
    cp0.ControlPointIndex = 0
    cp0.CumulativeMetersetWeight = 0.0
    cp0.GantryAngle = 0.0
    cp0.GantryRotationDirection = "NONE"
    cp0.BeamLimitingDeviceAngle = 0.0
    cp0.BeamLimitingDeviceRotationDirection = "NONE"
    cp0.PatientSupportAngle = 0.0
    cp0.PatientSupportRotationDirection = "NONE"
    cp0.NominalBeamEnergy = 6.0
    cp0.IsocenterPosition = [1.5, 1.5, 2.5]
    cp0.BeamLimitingDevicePositionSequence = Sequence([cp0_position])

    cp1 = Dataset()
    cp1.ControlPointIndex = 1
    cp1.CumulativeMetersetWeight = 1.0

    beam = Dataset()
    beam.BeamNumber = 1
    beam.BeamName = "B1"
    beam.BeamType = "STATIC"
    beam.RadiationType = "PHOTON"
    beam.TreatmentDeliveryType = "TREATMENT"
    beam.SourceAxisDistance = 1000.0
    beam.FinalCumulativeMetersetWeight = 1.0
    beam.NumberOfControlPoints = 2
    beam.NumberOfWedges = 0
    beam.NumberOfCompensators = 0
    beam.NumberOfBoli = 0
    beam.NumberOfBlocks = 0
    beam.BeamLimitingDeviceSequence = Sequence([device])
    beam.ControlPointSequence = Sequence([cp0, cp1])
    ds.BeamSequence = Sequence([beam])

    referenced_beam = Dataset()
    referenced_beam.ReferencedBeamNumber = 1
    referenced_beam.BeamMeterset = 100.0

    fraction = Dataset()
    fraction.FractionGroupNumber = 1
    fraction.NumberOfFractionsPlanned = 1
    fraction.ReferencedBeamSequence = Sequence([referenced_beam])
    ds.FractionGroupSequence = Sequence([fraction])

    ds.save_as(path, enforce_file_format=True)


def _build_case(root: Path) -> None:
    study_uid = generate_uid()
    series_uid = generate_uid()
    frame_uid = generate_uid()
    z_values = (0.0, 2.5, 5.0)

    sops = tuple(
        _write_ct_slice(
            root / f"ct_{index}.dcm",
            series_uid=series_uid,
            study_uid=study_uid,
            frame_uid=frame_uid,
            z_mm=z,
        )
        for index, z in enumerate(z_values)
    )

    struct_uid = _write_rtstruct(
        root / "struct.dcm",
        study_uid=study_uid,
        frame_uid=frame_uid,
        ct_series_uid=series_uid,
        ct_sops=sops,
        z_values=z_values,
    )
    _write_rtplan(
        root / "plan.dcm",
        study_uid=study_uid,
        struct_uid=struct_uid,
    )


def test_build_patient_artifacts_end_to_end(tmp_path: Path) -> None:
    dicom = tmp_path / "dicom"
    output = tmp_path / "out"
    dicom.mkdir()
    _build_case(dicom)

    result = build_patient_artifacts(
        dicom,
        output,
        patient_roi_name="patient",
        mixture_bins=4,
        out_of_range="raise",
    )

    assert result.egsphant_path.exists()
    assert result.media_definition_path.exists()
    assert result.summary_path.exists()
    assert result.medium_count >= 1
    assert result.low_hu_count == 0
    assert result.high_hu_count == 0

    phantom = read_egsphant(result.egsphant_path)
    assert phantom.geometry.shape_zyx == (3, 4, 4)
    assert phantom.medium_index.shape == (3, 4, 4)

    summary = json.loads(result.summary_path.read_text(encoding="utf-8"))
    assert summary["ct"]["patient_position"] == "HFS"
    assert summary["ct"]["kvp"] == 120.0
    assert summary["patient_roi"]["name"] == "patient"
    assert summary["patient_model"]["mixture_bins"] == 4
    assert summary["privacy"]["patient_name_stored"] is False
    assert summary["privacy"]["patient_id_stored"] is False
    assert len(summary["artifacts"]["egsphant"]["sha256"]) == 64


def test_patient_build_refuses_overwrite(tmp_path: Path) -> None:
    dicom = tmp_path / "dicom"
    output = tmp_path / "out"
    dicom.mkdir()
    _build_case(dicom)

    build_patient_artifacts(
        dicom,
        output,
        patient_roi_name="patient",
        mixture_bins=4,
    )

    try:
        build_patient_artifacts(
            dicom,
            output,
            patient_roi_name="patient",
            mixture_bins=4,
        )
    except FileExistsError:
        pass
    else:
        raise AssertionError("Existing patient artifacts were overwritten silently.")



def test_build_patient_artifacts_with_explicit_ct_derived_mask(
    tmp_path: Path,
) -> None:
    dicom = tmp_path / "dicom-derived"
    output = tmp_path / "out-derived"
    dicom.mkdir()
    _build_case(dicom)

    result = build_patient_artifacts(
        dicom,
        output,
        derive_external_from_isocenter=True,
        external_threshold_hu=-500.0,
        external_closing_iterations=0,
        mixture_bins=4,
        out_of_range="raise",
    )

    summary = json.loads(result.summary_path.read_text(encoding="utf-8"))
    assert summary["patient_mask"]["source_type"] == "CT_DERIVED_RESEARCH"
    assert summary["patient_mask"]["rtstruct_roi"] is None
    derived = summary["patient_mask"]["ct_derived"]
    assert derived["validation_status"] == "RESEARCH_ONLY"
    assert derived["seed_patient_mm"] == [1.5, 1.5, 2.5]
