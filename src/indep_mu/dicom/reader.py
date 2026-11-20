"""Conservative DICOM reader for standard Part 10 and raw Monaco datasets.

Some Monaco exports contain valid implicit-VR little-endian datasets without
the 128-byte preamble, DICM prefix or file-meta information. In these cases
strict pydicom.dcmread(force=False) rejects the file. This reader supports
that one legacy encoding only; it must not silently accept arbitrary bytes as
a DICOM object.
"""

from __future__ import annotations

from pathlib import Path
import struct

import pydicom
from pydicom.dataset import FileMetaDataset
from pydicom.errors import InvalidDicomError
from pydicom.uid import ImplicitVRLittleEndian


def _looks_like_raw_implicit_dataset(path: Path) -> bool:
    with path.open("rb") as stream:
        header = stream.read(16)
    if len(header) < 16:
        return False

    group, element, length = struct.unpack_from("<HHI", header, 0)
    # Native Monaco CT/RT exports begin with the implicit-VR study group;
    # do not activate force=True for arbitrary unsupported formats.
    if group != 0x0008 or element not in {0x0000, 0x0005, 0x0008, 0x0012, 0x0016}:
        return False
    return 0 <= length < 1_000_000


def dcmread_rt(
    path: str | Path,
    *,
    stop_before_pixels: bool = False,
) -> pydicom.dataset.FileDataset:
    """Read a Part 10 file or a validated headerless implicit-VR dataset."""

    source = Path(path)
    try:
        return pydicom.dcmread(
            source,
            stop_before_pixels=stop_before_pixels,
            force=False,
        )
    except InvalidDicomError:
        if not _looks_like_raw_implicit_dataset(source):
            raise ValueError(
                "DICOM file is not Part 10 and does not resemble the supported "
                "raw implicit-VR little-endian Monaco export."
            ) from None

    ds = pydicom.dcmread(
        source,
        stop_before_pixels=stop_before_pixels,
        force=True,
    )

    # Never interpret unknown binary data as patient DICOM merely because
    # force=True happens to parse a few tags.
    expected_classes = {
        "CT": "1.2.840.10008.5.1.4.1.1.2",
        "RTSTRUCT": "1.2.840.10008.5.1.4.1.1.481.3",
        "RTPLAN": "1.2.840.10008.5.1.4.1.1.481.5",
        "RTDOSE": "1.2.840.10008.5.1.4.1.1.481.2",
    }
    modality = str(getattr(ds, "Modality", ""))
    sop_class = str(getattr(ds, "SOPClassUID", ""))
    uid = str(getattr(ds, "SOPInstanceUID", ""))
    if modality not in expected_classes or sop_class != expected_classes[modality] or not uid:
        raise ValueError(
            "Raw DICOM must have an approved CT/RT modality, matching "
            "SOPClassUID and a non-empty SOPInstanceUID."
        )

    if not hasattr(ds, "file_meta") or ds.file_meta is None:
        ds.file_meta = FileMetaDataset()
    ds.file_meta.TransferSyntaxUID = ImplicitVRLittleEndian
    return ds
