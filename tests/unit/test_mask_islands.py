import numpy as np
from scipy import ndimage

from indep_mu.patient_model.mask_qc import survey_posterior_mask_islands


def _connected_body_and_support() -> np.ndarray:
    mask = np.zeros((5, 20, 40), dtype=bool)
    mask[:, 3:11, 10:30] = True
    mask[:, 15:17, 5:35] = True
    # One thin attachment makes all support and body voxels 3-D connected,
    # while the two regions are disconnected on four axial slices.
    mask[2, 10:16, 12] = True
    return mask


def test_posterior_islands_detect_included_support() -> None:
    mask = _connected_body_and_support()
    _, count3d = ndimage.label(
        mask,
        structure=ndimage.generate_binary_structure(3, 1),
    )
    assert count3d == 1

    report = survey_posterior_mask_islands(mask)
    assert report.suspicious_slice_indices == (0, 1, 3, 4)
    assert report.suspicious_voxel_count == 4 * 2 * 30


def test_no_posterior_islands_is_not_independent_clearance() -> None:
    mask = np.zeros((3, 20, 40), dtype=bool)
    mask[:, 3:11, 10:30] = True
    result = survey_posterior_mask_islands(mask)
    assert result.suspicious_voxel_count == 0
    assert result.suspicious_slice_indices == ()


def test_small_posterior_noise_is_suppressed() -> None:
    mask = np.zeros((3, 20, 40), dtype=bool)
    mask[:, 3:11, 10:30] = True
    mask[:, 16, 2:8] = True
    result = survey_posterior_mask_islands(mask)
    assert result.suspicious_voxel_count == 0
