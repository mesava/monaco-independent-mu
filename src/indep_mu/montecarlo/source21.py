from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol, runtime_checkable

import numpy as np

from indep_mu.dicom.rtplan import Beam
from indep_mu.patient_model.egsphant import EgsphantGeometry

from .syncjaws import SyncJawsSequence
from .syncmlce import SyncMlceSequence


@dataclass(frozen=True)
class Source21Angles:
    theta_deg: float
    phi_deg: float
    phicol_deg: float

    def __post_init__(self) -> None:
        values = np.asarray(
            [self.theta_deg, self.phi_deg, self.phicol_deg],
            dtype=np.float64,
        )
        if not np.all(np.isfinite(values)):
            raise ValueError("Source 21 angles must be finite.")


@runtime_checkable
class Source21OrientationMapper(Protocol):
    """Explicit IEC/DICOM -> DOSXYZnrc orientation transform.

    There is intentionally no default implementation here.  Gantry/couch/
    collimator conversion is non-trivial and must be independently verified
    before a clinical mapper is selected.
    """

    def map_angles(
        self,
        *,
        gantry_deg: float,
        collimator_deg: float,
        patient_support_deg: float,
    ) -> Source21Angles:
        ...


@dataclass(frozen=True)
class Source21ControlPoint:
    xiso_cm: float
    yiso_cm: float
    ziso_cm: float
    theta_deg: float
    phi_deg: float
    phicol_deg: float
    dsource_cm: float
    mu_index: float


@dataclass(frozen=True)
class Source21ControlPointSet:
    beam_number: int
    points: tuple[Source21ControlPoint, ...]


def patient_point_to_egsphant_cm(
    point_patient_mm: tuple[float, float, float],
    geometry: EgsphantGeometry,
) -> tuple[float, float, float]:
    """Express a DICOM patient point in the egsphant axis basis."""

    point = np.asarray(point_patient_mm, dtype=np.float64)
    if point.shape != (3,) or not np.all(np.isfinite(point)):
        raise ValueError("Patient point must contain three finite coordinates.")

    axes = (
        np.asarray(geometry.x_direction_patient, dtype=np.float64),
        np.asarray(geometry.y_direction_patient, dtype=np.float64),
        np.asarray(geometry.z_direction_patient, dtype=np.float64),
    )
    result = tuple(float(np.dot(point, axis) / 10.0) for axis in axes)
    return result


def build_source21_control_points(
    beam: Beam,
    *,
    phantom_geometry: EgsphantGeometry,
    orientation_mapper: Source21OrientationMapper,
    dsource_cm: float,
    equal_mu_tolerance: float = 1e-10,
) -> Source21ControlPointSet:
    """Build DOSXYZnrc source-21 control points for one beam.

    dsource_cm is explicit because it belongs to the BEAM/DOSXYZ transport
    geometry, not to DICOM RTPLAN alone.
    """

    if dsource_cm <= 0:
        raise ValueError("dsource_cm must be positive.")
    if beam.final_cumulative_meterset_weight <= 0:
        raise ValueError("FinalCumulativeMetersetWeight must be positive.")

    points: list[Source21ControlPoint] = []

    for cp in beam.control_points:
        if cp.isocenter_position_mm is None:
            raise ValueError(
                f"Beam {beam.number} CP {cp.index} has no IsocenterPosition."
            )
        if cp.gantry_angle_deg is None:
            raise ValueError(
                f"Beam {beam.number} CP {cp.index} has no GantryAngle."
            )
        if cp.collimator_angle_deg is None:
            raise ValueError(
                f"Beam {beam.number} CP {cp.index} has no collimator angle."
            )
        if cp.patient_support_angle_deg is None:
            raise ValueError(
                f"Beam {beam.number} CP {cp.index} has no PatientSupportAngle."
            )

        xiso, yiso, ziso = patient_point_to_egsphant_cm(
            cp.isocenter_position_mm,
            phantom_geometry,
        )
        angles = orientation_mapper.map_angles(
            gantry_deg=cp.gantry_angle_deg,
            collimator_deg=cp.collimator_angle_deg,
            patient_support_deg=cp.patient_support_angle_deg,
        )
        point = Source21ControlPoint(
            xiso_cm=xiso,
            yiso_cm=yiso,
            ziso_cm=ziso,
            theta_deg=float(angles.theta_deg),
            phi_deg=float(angles.phi_deg),
            phicol_deg=float(angles.phicol_deg),
            dsource_cm=float(dsource_cm),
            mu_index=(
                cp.cumulative_meterset_weight
                / beam.final_cumulative_meterset_weight
            ),
        )

        if points and abs(point.mu_index - points[-1].mu_index) <= equal_mu_tolerance:
            points[-1] = point
        else:
            points.append(point)

    indices = np.asarray([point.mu_index for point in points], dtype=np.float64)
    if len(points) < 2:
        raise ValueError("Source 21 requires at least two control points.")
    if abs(indices[0]) > equal_mu_tolerance:
        raise ValueError("Source 21 control points must start at muIndex=0.")
    if abs(indices[-1] - 1.0) > equal_mu_tolerance:
        raise ValueError("Source 21 control points must end at muIndex=1.")
    if np.any(np.diff(indices) <= equal_mu_tolerance):
        raise ValueError("Source 21 muIndex values must be strictly increasing.")

    return Source21ControlPointSet(
        beam_number=beam.number,
        points=tuple(points),
    )


def render_source21_control_points(control_points: Source21ControlPointSet) -> str:
    """Render only SC1-21a records, one 8-value line per control point."""

    lines: list[str] = []
    for point in control_points.points:
        values = (
            point.xiso_cm,
            point.yiso_cm,
            point.ziso_cm,
            point.theta_deg,
            point.phi_deg,
            point.phicol_deg,
            point.dsource_cm,
            point.mu_index,
        )
        lines.append("".join(f"{value:15.8f}" for value in values))
    return "\n".join(lines) + "\n"


def validate_synchronized_mu_indices(
    source21: Source21ControlPointSet,
    *,
    mlc: SyncMlceSequence | None = None,
    jaws: SyncJawsSequence | None = None,
    tolerance: float = 1e-8,
) -> None:
    """Require source motion and synchronized component motion to share MUINDEX."""

    reference = np.asarray(
        [point.mu_index for point in source21.points],
        dtype=np.float64,
    )

    if mlc is not None:
        candidate = np.asarray(
            [point.mu_index for point in mlc.points],
            dtype=np.float64,
        )
        if candidate.shape != reference.shape or not np.allclose(
            candidate,
            reference,
            rtol=0.0,
            atol=tolerance,
        ):
            raise ValueError(
                "Source 21 and SYNCMLCE MUINDEX grids are not synchronized."
            )

    if jaws is not None:
        candidate = np.asarray(
            [point.mu_index for point in jaws.points],
            dtype=np.float64,
        )
        if candidate.shape != reference.shape or not np.allclose(
            candidate,
            reference,
            rtol=0.0,
            atol=tolerance,
        ):
            raise ValueError(
                "Source 21 and SYNCJAWS MUINDEX grids are not synchronized."
            )
