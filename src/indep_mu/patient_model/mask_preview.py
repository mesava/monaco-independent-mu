"""Local visual review of CT-derived external masks.

This module produces a private PNG for a qualified reviewer. It must not be
used to automatically label/clear support or to commit patient CT pixels.
Matplotlib is an optional dependency and is imported only when rendering.
"""

from __future__ import annotations

from pathlib import Path
import numpy as np

from indep_mu.dicom.ct import CtSeries


def render_external_mask_montage(
    ct: CtSeries,
    mask: np.ndarray,
    output_png: str | Path,
    *,
    panel_count: int = 9,
    priority_slice_indices: tuple[int, ...] = (),
    window_min_hu: float = -1000.0,
    window_max_hu: float = 700.0,
    overwrite: bool = False,
) -> tuple[int, ...]:
    """Render selected full-axial slices with a mask boundary in a local PNG.

    Returns the original CT z indices chosen for display. The rendered image
    contains no patient identifiers but remains patient-derived medical data.
    """

    if panel_count < 1 or panel_count > 20:
        raise ValueError("panel_count must be in [1, 20].")
    if not (np.isfinite(window_min_hu) and np.isfinite(window_max_hu)):
        raise ValueError("CT window limits must be finite.")
    if window_max_hu <= window_min_hu:
        raise ValueError("window_max_hu must exceed window_min_hu.")

    selected = np.asarray(mask, dtype=bool)
    if selected.shape != ct.hu.shape or selected.ndim != 3:
        raise ValueError("Mask shape must match the CT (z,y,x) image stack.")

    active = np.flatnonzero(np.any(selected, axis=(1, 2)))
    if active.size == 0:
        raise ValueError("Cannot render an empty patient mask.")
    output = Path(output_png)
    if output.exists() and not overwrite:
        raise FileExistsError(f"Output image already exists: {output}")

    # Preserve broad anatomical coverage, but also prioritize slices with
    # posterior 2-D islands detected inside the connected 3-D mask.
    # The islands are not automatically classified as couch.
    baseline_count = panel_count if not priority_slice_indices else max(
        1, panel_count // 2
    )
    baseline_sample = np.linspace(
        0, active.size - 1, min(baseline_count, active.size)
    )
    chosen = {
        int(active[int(np.rint(value))])
        for value in baseline_sample
    }

    if priority_slice_indices:
        candidates = sorted(
            set(int(value) for value in priority_slice_indices)
            & set(int(value) for value in active)
        )
        remaining = max(0, panel_count - len(chosen))
        if candidates and remaining:
            sample = np.linspace(0, len(candidates) - 1, min(remaining, len(candidates)))
            chosen.update(candidates[int(np.rint(value))] for value in sample)

    if len(chosen) < min(panel_count, active.size):
        candidates = [int(value) for value in active if int(value) not in chosen]
        remaining = min(panel_count, active.size) - len(chosen)
        if candidates and remaining > 0:
            sample = np.linspace(0, len(candidates) - 1, remaining)
            chosen.update(candidates[int(np.rint(value))] for value in sample)

    indices = tuple(sorted(chosen))

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise RuntimeError(
            "Local QC PNG requires matplotlib: pip install -e '.[qc]'."
        ) from exc

    columns = 3 if len(indices) > 3 else len(indices)
    rows = int(np.ceil(len(indices) / columns))
    fig, axes = plt.subplots(
        rows,
        columns,
        figsize=(5.0 * columns, 5.2 * rows),
        squeeze=False,
    )

    for ax in axes.flat:
        ax.set_visible(False)

    for ax, z_index in zip(axes.flat, indices):
        ax.set_visible(True)
        ax.imshow(
            ct.hu[z_index],
            cmap="gray",
            vmin=window_min_hu,
            vmax=window_max_hu,
            interpolation="nearest",
        )
        ax.contour(
            selected[z_index].astype(np.uint8),
            levels=(0.5,),
            colors=("red",),
            linewidths=0.8,
        )
        ax.set_title(
            f"Slice {z_index} | position {ct.slice_positions_mm[z_index]:.1f} mm"
        )
        ax.set_xlabel("DICOM column")
        ax.set_ylabel("DICOM row (posterior for axial HFS)")
        ax.set_xlim(0, ct.hu.shape[2] - 1)
        ax.set_ylim(ct.hu.shape[1] - 1, 0)

    fig.suptitle(
        "RESEARCH ONLY — CT-derived mask / support NOT independently cleared"
    )
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    try:
        fig.savefig(output, dpi=130)
    finally:
        plt.close(fig)

    return indices
