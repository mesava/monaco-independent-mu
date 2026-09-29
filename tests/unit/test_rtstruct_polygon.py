from __future__ import annotations

import numpy as np

from indep_mu.dicom.rtstruct import _polygon_mask


def test_polygon_mask_fills_square_at_voxel_centres() -> None:
    mask = _polygon_mask(
        8,
        8,
        row_coordinates=np.asarray([2.0, 2.0, 5.0, 5.0]),
        column_coordinates=np.asarray([2.0, 5.0, 5.0, 2.0]),
    )

    assert mask[3, 3]
    assert mask[4, 4]
    assert not mask[1, 3]
    assert not mask[6, 3]


def test_polygon_mask_handles_concave_polygon() -> None:
    mask = _polygon_mask(
        8,
        8,
        row_coordinates=np.asarray([1, 1, 6, 6, 3, 3]),
        column_coordinates=np.asarray([1, 6, 6, 3, 3, 1]),
    )

    assert mask[2, 5]
    assert mask[5, 4]
    assert not mask[5, 2]
