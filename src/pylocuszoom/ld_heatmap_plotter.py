"""LD heatmap generator for pairwise linkage disequilibrium visualization.

Provides triangular heatmap display of pairwise LD values (R² or D')
with colorbar legend and SNP highlighting support.
"""

from typing import Any, List, Optional, Tuple, Union

import numpy as np
import pandas as pd

from ._figure import FigurePlan, render_figure
from .backends import BackendType, get_backend
from .config import LDMetric
from .panels.ld_heatmap import LDHeatmapPanel


class LDHeatmapPlotter:
    """LD heatmap generator for pairwise LD visualization.

    Creates triangular heatmaps showing pairwise linkage disequilibrium
    between variants. Supports R² and D' metrics, lead SNP highlighting,
    and multiple backend renderers.

    Supports multiple rendering backends:
    - matplotlib (default): Static publication-quality plots
    - plotly: Interactive HTML with hover tooltips
    - bokeh: Interactive HTML for dashboards

    Args:
        backend: Plotting backend ('matplotlib', 'plotly', or 'bokeh').

    Example:
        >>> plotter = LDHeatmapPlotter()
        >>> fig = plotter.plot_ld_heatmap(ld_matrix, lead_snp="rs12345")
        >>> fig.savefig("ld_heatmap.png", dpi=150)
    """

    def __init__(
        self,
        backend: BackendType = "matplotlib",
    ):
        """Initialize the LD heatmap plotter."""
        self._backend = get_backend(backend)

    def plot_ld_heatmap(
        self,
        ld_matrix: Union[pd.DataFrame, np.ndarray],
        snp_ids: Optional[List[str]] = None,
        lead_snp: Optional[str] = None,
        highlight_snps: Optional[List[str]] = None,
        metric: LDMetric = "r2",
        figsize: Tuple[float, float] = (8, 8),
        title: Optional[str] = None,
        show_colorbar: bool = True,
    ) -> Any:
        """Create triangular LD heatmap.

        Args:
            ld_matrix: Square DataFrame or numpy array with pairwise LD values.
                NaN values are displayed as grey (missing data).
            snp_ids: List of SNP IDs for axis labels. If None, uses matrix index.
            lead_snp: SNP ID to highlight as lead variant (red highlight).
            highlight_snps: Additional SNP IDs to highlight (blue highlight).
            metric: LD metric label for colorbar ("r2" or "dprime").
            figsize: Figure size as (width, height).
            title: Plot title.
            show_colorbar: Whether to show colorbar legend.

        Returns:
            Figure object (type depends on backend).

        Raises:
            ValidationError: If ld_matrix is not square, lead_snp or any
                highlight_snps is not in snp_ids, or metric is not "r2" or
                "dprime".

        Example:
            >>> fig = plotter.plot_ld_heatmap(
            ...     ld_matrix,
            ...     snp_ids=["rs1", "rs2", "rs3"],
            ...     lead_snp="rs1",
            ...     metric="r2",
            ... )
        """
        panel = LDHeatmapPanel.from_matrix(
            ld_matrix,
            snp_ids,
            lead_snp=lead_snp,
            highlight_snps=highlight_snps,
            metric=metric,
            title=title,
            show_colorbar=show_colorbar,
        )
        return render_figure(self._backend, FigurePlan(panels=[panel], figsize=figsize))
