"""LD heatmap generator for pairwise linkage disequilibrium visualization.

Provides triangular heatmap display of pairwise LD values (R² or D')
with colorbar legend and SNP highlighting support.
"""

from typing import Any, List, Optional, Tuple, Union

import numpy as np
import pandas as pd

from ._ld_heatmap_renderer import LDHeatmapRequest, render_ld_heatmap
from .backends import BackendType, get_backend


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
        metric: str = "r2",
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
            ValueError: If ld_matrix is not square.
            ValueError: If lead_snp not found in snp_ids.
            ValueError: If any highlight_snps not found in snp_ids.

        Example:
            >>> fig = plotter.plot_ld_heatmap(
            ...     ld_matrix,
            ...     snp_ids=["rs1", "rs2", "rs3"],
            ...     lead_snp="rs1",
            ...     metric="r2",
            ... )
        """
        # Extract data and snp_ids from DataFrame if needed
        if isinstance(ld_matrix, pd.DataFrame):
            data = ld_matrix.values
            if snp_ids is None:
                snp_ids = list(ld_matrix.index.astype(str))
        else:
            data = np.asarray(ld_matrix)
            if snp_ids is None:
                snp_ids = [str(i) for i in range(data.shape[0])]

        # Validate square matrix
        if data.ndim != 2 or data.shape[0] != data.shape[1]:
            raise ValueError(f"ld_matrix must be square, got shape {data.shape}")

        n_snps = len(snp_ids)
        if data.shape[0] != n_snps:
            raise ValueError(
                f"snp_ids length ({n_snps}) does not match matrix dimension ({data.shape[0]})"
            )

        # Validate lead_snp
        lead_idx = None
        if lead_snp is not None:
            if lead_snp not in snp_ids:
                raise ValueError(f"lead_snp '{lead_snp}' not found in snp_ids")
            lead_idx = snp_ids.index(lead_snp)

        # Validate highlight_snps
        highlight_indices = []
        if highlight_snps:
            for snp in highlight_snps:
                if snp not in snp_ids:
                    raise ValueError(f"highlight_snp '{snp}' not found in snp_ids")
                highlight_indices.append(snp_ids.index(snp))

        return render_ld_heatmap(
            self._backend,
            LDHeatmapRequest(
                data=data,
                snp_ids=snp_ids,
                lead_idx=lead_idx,
                highlight_indices=highlight_indices,
                metric=metric,
                figsize=figsize,
                title=title,
                show_colorbar=show_colorbar,
            ),
        )
