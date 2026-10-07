from __future__ import annotations

from pathlib import Path

import pytest
from pydicom.dataset import Dataset, FileDataset, FileMetaDataset
from pydicom.sequence import Sequence
from pydicom.uid import ExplicitVRLittleEndian, RTPlanStorage, generate_uid

from indep_mu.dicom.rtplan import load_rtplan


def _device_definition(
    device_type: str,
    pairs: int,
    boundaries: list[float] | None = None,
) -> Dataset:
    item = Dataset()
    item.RTBeamLimitingDeviceType = device_type
    item.SourceToBeamLimitingDeviceDistance = 500.0
    item.NumberOfLeafJawPairs = pairs
    if boundaries is not None:
        item.LeafPositionBoundaries = boundaries
    return item


def _device_state(device_type: str, positions: list[float]) -> Dataset:
    item = Dataset()
    item.RTBeamLimitingDeviceType = device_type
    item.LeafJawPositions = positions
    return item


def _control_point(
    index: int,
    cmw: float,
    *,
    gantry: float | None = None,
    jaws: list[float] | None = None,
    mlc: list[float] | None = None,
) -> Dataset:
    cp = Dataset()
    cp.ControlPointIndex = index
    cp.CumulativeMetersetWeight = cmw

    if gantry is not None:
        cp.GantryAngle = gantry
        cp.GantryRotationDirection = "CC"

    states = []
    if jaws is not None:
        states.append(_device_state("ASYMY", jaws))
    if mlc is not None:
        states.append(_device_state("MLCX", mlc))
    if states:
        cp.BeamLimitingDevicePositionSequence = Sequence(states)

    if index == 0:
        cp.BeamLimitingDeviceAngle = 15.0
        cp.BeamLimitingDeviceRotationDirection = "NONE"
        cp.PatientSupportAngle = 0.0
        cp.PatientSupportRotationDirection = "NONE"
        cp.IsocenterPosition = [10.0, 20.0, 30.0]
        cp.NominalBeamEnergy = 6.0
        cp.DoseRateSet = 600.0
        cp.SourceToSurfaceDistance = 900.0

    return cp


def _write_plan(path: Path, cmw_values=(0.0, 0.4, 1.0)) -> None:
    meta = FileMetaDataset()
    meta.MediaStorageSOPClassUID = RTPlanStorage
    meta.MediaStorageSOPInstanceUID = generate_uid()
    meta.TransferSyntaxUID = ExplicitVRLittleEndian

    ds = FileDataset(str(path), {}, file_meta=meta, preamble=b"\0" * 128)
    ds.SOPClassUID = RTPlanStorage
    ds.SOPInstanceUID = meta.MediaStorageSOPInstanceUID
    ds.Modality = "RTPLAN"
    ds.RTPlanLabel = "TEST"
    ds.RTPlanName = "Synthetic plan"

    fraction = Dataset()
    fraction.FractionGroupNumber = 1
    fraction.NumberOfFractionsPlanned = 15
    referenced = Dataset()
    referenced.ReferencedBeamNumber = 1
    referenced.BeamMeterset = 500.0
    fraction.ReferencedBeamSequence = Sequence([referenced])
    ds.FractionGroupSequence = Sequence([fraction])

    beam = Dataset()
    beam.BeamNumber = 1
    beam.BeamName = "ARC1"
    beam.BeamDescription = "Synthetic VMAT"
    beam.BeamType = "DYNAMIC"
    beam.RadiationType = "PHOTON"
    beam.TreatmentDeliveryType = "TREATMENT"
    beam.TreatmentMachineName = "VERSA"
    beam.SourceAxisDistance = 1000.0
    beam.FinalCumulativeMetersetWeight = 1.0

    fluence = Dataset()
    fluence.FluenceMode = "STANDARD"
    beam.PrimaryFluenceModeSequence = Sequence([fluence])

    beam.BeamLimitingDeviceSequence = Sequence(
        [
            _device_definition("ASYMY", 1),
            _device_definition("MLCX", 2, [-10.0, -5.0, 0.0]),
        ]
    )

    cps = [
        _control_point(
            0,
            cmw_values[0],
            gantry=180.0,
            jaws=[-50.0, 50.0],
            mlc=[-10.0, -8.0, 8.0, 10.0],
        ),
        _control_point(
            1,
            cmw_values[1],
            gantry=150.0,
            mlc=[-9.0, -7.0, 7.0, 9.0],
        ),
        _control_point(
            2,
            cmw_values[2],
            gantry=120.0,
            mlc=[-8.0, -6.0, 6.0, 8.0],
        ),
    ]
    beam.NumberOfControlPoints = len(cps)
    beam.ControlPointSequence = Sequence(cps)
    ds.BeamSequence = Sequence([beam])

    ds.save_as(path, enforce_file_format=True)


def test_rtplan_parses_meterset_segments_and_inheritance(tmp_path: Path) -> None:
    path = tmp_path / "plan.dcm"
    _write_plan(path)

    plan = load_rtplan(path)
    beam = plan.beam(1)

    assert plan.number_of_fractions_planned == 15
    assert beam.beam_meterset_mu == pytest.approx(500.0)
    assert beam.final_cumulative_meterset_weight == pytest.approx(1.0)
    assert beam.control_points[1].positions("ASYMY") == (-50.0, 50.0)
    assert beam.control_points[1].nominal_beam_energy_mv == pytest.approx(6.0)
    assert beam.control_points[2].isocenter_position_mm == (10.0, 20.0, 30.0)

    segments = beam.segments
    assert len(segments) == 2
    assert segments[0].delta_mu == pytest.approx(200.0)
    assert segments[1].delta_mu == pytest.approx(300.0)
    assert sum(segment.delta_mu for segment in segments) == pytest.approx(500.0)


def test_rtplan_rejects_nonmonotonic_cmw(tmp_path: Path) -> None:
    path = tmp_path / "plan.dcm"
    _write_plan(path, cmw_values=(0.0, 0.6, 0.5))

    with pytest.raises(ValueError, match="not monotonic|FinalCumulative"):
        load_rtplan(path)


def test_rtplan_rejects_missing_initial_device_state(tmp_path: Path) -> None:
    path = tmp_path / "plan.dcm"
    _write_plan(path)

    import pydicom

    ds = pydicom.dcmread(path)
    first = ds.BeamSequence[0].ControlPointSequence[0]
    first.BeamLimitingDevicePositionSequence = Sequence(
        [first.BeamLimitingDevicePositionSequence[1]]
    )
    ds.save_as(path, enforce_file_format=True)

    with pytest.raises(ValueError, match="no inherited state"):
        load_rtplan(path)
