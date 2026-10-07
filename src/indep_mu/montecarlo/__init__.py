"""Monte Carlo backend and dynamic-delivery abstractions."""

from .absolute_calibration import AbsoluteMcCalibration, assert_reference_geometry_matches
from .delivery import (
    DeliveryInterpolationPolicy,
    SampledBeamState,
    directed_rotation_delta_deg,
    sample_beam_segment,
)
from .interface import (
    MonteCarloBackend,
    MonteCarloCapabilities,
    MonteCarloRunConfig,
)
from .source21 import (
    Source21Angles,
    Source21ControlPoint,
    Source21ControlPointSet,
    Source21OrientationMapper,
    build_source21_control_points,
    patient_point_to_egsphant_cm,
    render_source21_control_points,
    validate_synchronized_mu_indices,
)
from .source21_reference import (
    CoplanarHfsReferenceAngles,
    coplanar_hfs_reference_angles,
)
from .syncjaws import (
    FocusedJawPairGeometry,
    SyncJawsPairState,
    SyncJawsSequence,
    SyncJawsSequencePoint,
    build_syncjaws_sequence,
    render_syncjaws_sequence,
    write_syncjaws_sequence,
)
from .syncmlce import (
    SourceFocusedSyncMlceMapper,
    SyncMlceCoordinateMapper,
    SyncMlceOpening,
    SyncMlceSequence,
    SyncMlceSequencePoint,
    build_syncmlce_sequence,
    render_syncmlce_sequence,
    write_syncmlce_sequence,
)
from .threedose import ThreeDDose, read_3ddose, write_3ddose

__all__ = [
    "AbsoluteMcCalibration",
    "CoplanarHfsReferenceAngles",
    "DeliveryInterpolationPolicy",
    "FocusedJawPairGeometry",
    "MonteCarloBackend",
    "MonteCarloCapabilities",
    "MonteCarloRunConfig",
    "SampledBeamState",
    "Source21Angles",
    "Source21ControlPoint",
    "Source21ControlPointSet",
    "Source21OrientationMapper",
    "SourceFocusedSyncMlceMapper",
    "SyncJawsPairState",
    "SyncJawsSequence",
    "SyncJawsSequencePoint",
    "SyncMlceCoordinateMapper",
    "SyncMlceOpening",
    "SyncMlceSequence",
    "SyncMlceSequencePoint",
    "ThreeDDose",
    "assert_reference_geometry_matches",
    "build_source21_control_points",
    "build_syncjaws_sequence",
    "build_syncmlce_sequence",
    "coplanar_hfs_reference_angles",
    "directed_rotation_delta_deg",
    "patient_point_to_egsphant_cm",
    "read_3ddose",
    "render_source21_control_points",
    "render_syncjaws_sequence",
    "render_syncmlce_sequence",
    "sample_beam_segment",
    "validate_synchronized_mu_indices",
    "write_3ddose",
    "write_syncjaws_sequence",
    "write_syncmlce_sequence",
]
