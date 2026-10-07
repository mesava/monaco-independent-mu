from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pydicom


@dataclass(frozen=True)
class DicomObjectRef:
    path: Path
    modality: str
    sop_instance_uid: str
    series_instance_uid: str | None
    frame_of_reference_uid: str | None


@dataclass(frozen=True)
class CtSeriesRef:
    series_instance_uid: str
    frame_of_reference_uid: str
    instances: tuple[DicomObjectRef, ...]


@dataclass(frozen=True)
class CaseMessage:
    severity: str
    code: str
    message: str


@dataclass(frozen=True)
class DicomCaseManifest:
    root: Path
    ct_series: CtSeriesRef
    rtstruct: DicomObjectRef
    rtplan: DicomObjectRef
    rtdoses: tuple[DicomObjectRef, ...]
    messages: tuple[CaseMessage, ...]

    @property
    def errors(self) -> tuple[CaseMessage, ...]:
        return tuple(item for item in self.messages if item.severity == "ERROR")

    @property
    def warnings(self) -> tuple[CaseMessage, ...]:
        return tuple(item for item in self.messages if item.severity == "WARNING")

    @property
    def passed(self) -> bool:
        return not self.errors


def _frame_uid(dataset: pydicom.dataset.Dataset) -> str | None:
    value = getattr(dataset, "FrameOfReferenceUID", None)
    if value is not None:
        return str(value)

    frames: set[str] = set()
    for item in getattr(dataset, "ReferencedFrameOfReferenceSequence", []) or []:
        if hasattr(item, "FrameOfReferenceUID"):
            frames.add(str(item.FrameOfReferenceUID))

    for item in getattr(dataset, "StructureSetROISequence", []) or []:
        if hasattr(item, "ReferencedFrameOfReferenceUID"):
            frames.add(str(item.ReferencedFrameOfReferenceUID))

    if len(frames) == 1:
        return next(iter(frames))
    return None


def _scan_dicom_objects(root: Path) -> tuple[DicomObjectRef, ...]:
    result: list[DicomObjectRef] = []

    for path in sorted(item for item in root.rglob("*") if item.is_file()):
        try:
            ds = pydicom.dcmread(path, stop_before_pixels=True, force=False)
        except Exception:
            continue

        modality = str(getattr(ds, "Modality", ""))
        if modality not in {"CT", "RTSTRUCT", "RTPLAN", "RTDOSE"}:
            continue
        if not hasattr(ds, "SOPInstanceUID"):
            raise ValueError(f"{path}: DICOM object has no SOPInstanceUID.")

        result.append(
            DicomObjectRef(
                path=path,
                modality=modality,
                sop_instance_uid=str(ds.SOPInstanceUID),
                series_instance_uid=(
                    str(ds.SeriesInstanceUID)
                    if hasattr(ds, "SeriesInstanceUID")
                    else None
                ),
                frame_of_reference_uid=_frame_uid(ds),
            )
        )

    if not result:
        raise ValueError(f"No CT/RT DICOM objects found under {root}.")

    return tuple(result)


def _referenced_ct_series_uids(rtstruct_path: Path) -> set[str]:
    ds = pydicom.dcmread(rtstruct_path, stop_before_pixels=True, force=False)
    result: set[str] = set()

    for frame in getattr(ds, "ReferencedFrameOfReferenceSequence", []) or []:
        for study in getattr(frame, "RTReferencedStudySequence", []) or []:
            for series in getattr(study, "RTReferencedSeriesSequence", []) or []:
                if hasattr(series, "SeriesInstanceUID"):
                    result.add(str(series.SeriesInstanceUID))

    return result


def _referenced_ct_sop_uids(rtstruct_path: Path) -> set[str]:
    ds = pydicom.dcmread(rtstruct_path, stop_before_pixels=True, force=False)
    result: set[str] = set()

    for roi_contour in getattr(ds, "ROIContourSequence", []) or []:
        for contour in getattr(roi_contour, "ContourSequence", []) or []:
            for image in getattr(contour, "ContourImageSequence", []) or []:
                uid = getattr(image, "ReferencedSOPInstanceUID", None)
                if uid is not None:
                    result.add(str(uid))

    return result


def _rtplan_struct_uid(rtplan_path: Path) -> str | None:
    ds = pydicom.dcmread(rtplan_path, stop_before_pixels=True, force=False)
    sequence = getattr(ds, "ReferencedStructureSetSequence", None)
    if not sequence:
        return None
    if len(sequence) != 1:
        raise ValueError("RTPLAN references more than one RTSTRUCT.")
    return str(sequence[0].ReferencedSOPInstanceUID)


def _rtdose_plan_uids(rtdose_path: Path) -> set[str]:
    ds = pydicom.dcmread(rtdose_path, stop_before_pixels=True, force=False)
    return {
        str(item.ReferencedSOPInstanceUID)
        for item in getattr(ds, "ReferencedRTPlanSequence", []) or []
        if hasattr(item, "ReferencedSOPInstanceUID")
    }


def discover_dicom_case(path: str | Path) -> DicomCaseManifest:
    """Discover one strict CT + RTSTRUCT + RTPLAN case without patient identity."""

    root = Path(path)
    if not root.exists() or not root.is_dir():
        raise ValueError("DICOM case path must be an existing directory.")

    objects = _scan_dicom_objects(root)
    ct_objects = tuple(item for item in objects if item.modality == "CT")
    structs = tuple(item for item in objects if item.modality == "RTSTRUCT")
    plans = tuple(item for item in objects if item.modality == "RTPLAN")
    doses = tuple(item for item in objects if item.modality == "RTDOSE")

    if len(structs) != 1:
        raise ValueError(f"Expected exactly one RTSTRUCT, found {len(structs)}.")
    if len(plans) != 1:
        raise ValueError(f"Expected exactly one RTPLAN, found {len(plans)}.")
    if not ct_objects:
        raise ValueError("No CT instances found.")

    rtstruct = structs[0]
    rtplan = plans[0]
    messages: list[CaseMessage] = []

    by_series: dict[str, list[DicomObjectRef]] = {}
    for item in ct_objects:
        if item.series_instance_uid is None:
            raise ValueError(f"CT object {item.path.name} has no SeriesInstanceUID.")
        by_series.setdefault(item.series_instance_uid, []).append(item)

    referenced_series = _referenced_ct_series_uids(rtstruct.path)
    candidate_series = set(by_series) & referenced_series

    if len(candidate_series) == 1:
        selected_uid = next(iter(candidate_series))
    else:
        # Fallback remains reference-based: resolve contour-referenced CT SOPs.
        referenced_sops = _referenced_ct_sop_uids(rtstruct.path)
        matched = {
            series_uid
            for series_uid, instances in by_series.items()
            if referenced_sops
            and referenced_sops.issubset(
                {item.sop_instance_uid for item in instances}
            )
        }
        if len(matched) != 1:
            raise ValueError(
                "Cannot prove which CT series is referenced by RTSTRUCT."
            )
        selected_uid = next(iter(matched))

    selected_instances = tuple(by_series[selected_uid])
    frame_uids = {
        item.frame_of_reference_uid
        for item in selected_instances
        if item.frame_of_reference_uid is not None
    }
    if len(frame_uids) != 1:
        raise ValueError(
            "Selected CT series does not have exactly one FrameOfReferenceUID."
        )
    ct_frame_uid = next(iter(frame_uids))

    if rtstruct.frame_of_reference_uid is not None:
        if rtstruct.frame_of_reference_uid != ct_frame_uid:
            messages.append(
                CaseMessage(
                    "ERROR",
                    "RTSTRUCT_FRAME_MISMATCH",
                    "RTSTRUCT and referenced CT use different FrameOfReferenceUID.",
                )
            )
    else:
        messages.append(
            CaseMessage(
                "WARNING",
                "RTSTRUCT_FRAME_UNRESOLVED",
                "RTSTRUCT FrameOfReferenceUID could not be resolved uniquely.",
            )
        )

    struct_uid = _rtplan_struct_uid(rtplan.path)
    if struct_uid is None:
        messages.append(
            CaseMessage(
                "WARNING",
                "RTPLAN_NO_STRUCT_REFERENCE",
                "RTPLAN has no ReferencedStructureSetSequence.",
            )
        )
    elif struct_uid != rtstruct.sop_instance_uid:
        messages.append(
            CaseMessage(
                "ERROR",
                "RTPLAN_STRUCT_MISMATCH",
                "RTPLAN references a different RTSTRUCT SOP Instance UID.",
            )
        )

    for dose in doses:
        if (
            dose.frame_of_reference_uid is not None
            and dose.frame_of_reference_uid != ct_frame_uid
        ):
            messages.append(
                CaseMessage(
                    "ERROR",
                    "RTDOSE_FRAME_MISMATCH",
                    f"RTDOSE {dose.path.name} uses a different FrameOfReferenceUID.",
                )
            )

        plan_uids = _rtdose_plan_uids(dose.path)
        if plan_uids and rtplan.sop_instance_uid not in plan_uids:
            messages.append(
                CaseMessage(
                    "ERROR",
                    "RTDOSE_PLAN_MISMATCH",
                    f"RTDOSE {dose.path.name} references a different RTPLAN.",
                )
            )
        elif not plan_uids:
            messages.append(
                CaseMessage(
                    "WARNING",
                    "RTDOSE_PLAN_REFERENCE_MISSING",
                    f"RTDOSE {dose.path.name} has no ReferencedRTPlanSequence.",
                )
            )

    if not doses:
        messages.append(
            CaseMessage(
                "WARNING",
                "NO_RTDOSE",
                "No RTDOSE objects were found; TPS dose comparison is unavailable.",
            )
        )

    ct_series = CtSeriesRef(
        series_instance_uid=selected_uid,
        frame_of_reference_uid=ct_frame_uid,
        instances=selected_instances,
    )

    return DicomCaseManifest(
        root=root,
        ct_series=ct_series,
        rtstruct=rtstruct,
        rtplan=rtplan,
        rtdoses=doses,
        messages=tuple(messages),
    )
