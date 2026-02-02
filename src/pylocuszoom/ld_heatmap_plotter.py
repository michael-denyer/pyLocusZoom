"""LD heatmap generator for pairwise linkage disequilibrium visualization.

Provides triangular heatmap display of pairwise LD values (R² or D')
with colorbar legend and SNP highlighting support.
"""

from typing import Any, List, Optional, Tuple, Union

import numpy as np
import pandas as pd

from .backends import BackendType, get_backend
from .colors import (
    LD_HEATMAP_COLORS,
    LEAD_SNP_HIGHLIGHT_COLOR,
    SECONDARY_HIGHLIGHT_COLOR,
)


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
        species: Species name ('canine', 'feline', 'human', or None).
            Currently unused but kept for API consistency.
        backend: Plotting backend ('matplotlib', 'plotly', or 'bokeh').

    Example:
        >>> plotter = LDHeatmapPlotter()
        >>> fig = plotter.plot_ld_heatmap(ld_matrix, lead_snp="rs12345")
        >>> fig.savefig("ld_heatmap.png", dpi=150)
    """

    def __init__(
        self,
        species: str = "canine",
        backend: BackendType = "matplotlib",
    ):
        """Initialize the LD heatmap plotter."""
        self.species = species  # Kept for backward compatibility, currently unused
        self._backend = get_backend(backend)
        self.backend_name = backend

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

        # Create figure with single panel
        fig, axes = self._backend.create_figure(
            n_panels=1,
            height_ratios=[1.0],
            figsize=figsize,
            sharex=False,
        )
        ax = axes[0]

        # Render triangular heatmap
        mappable = self._backend.add_heatmap(
            ax,
            data=data,
            x_coords=list(range(n_snps)),
            y_coords=list(range(n_snps)),
            cmap_colors=LD_HEATMAP_COLORS,
            vmin=0.0,
            vmax=1.0,
            mask_upper=True,
        )

        # Add colorbar
        if show_colorbar:
            label = "R²" if metric == "r2" else "D'"
            self._backend.add_colorbar(ax, mappable, label=label)

        # Highlight lead SNP
        if lead_idx is not None:
            self._highlight_snp(
                ax=ax,
                fig=fig,
                snp_idx=lead_idx,
                n_snps=n_snps,
                color=LEAD_SNP_HIGHLIGHT_COLOR,
            )

        # Highlight additional SNPs
        for idx in highlight_indices:
            self._highlight_snp(
                ax=ax,
                fig=fig,
                snp_idx=idx,
                n_snps=n_snps,
                color=SECONDARY_HIGHLIGHT_COLOR,
            )

        # Set axis ticks with SNP labels
        tick_positions = list(range(n_snps))
        self._backend.set_xticks(ax, tick_positions, snp_ids, rotation=90)
        self._backend.set_yticks(ax, tick_positions, snp_ids)

        # Set title
        if title:
            self._backend.set_title(ax, title)

        # Finalize layout
        self._backend.finalize_layout(fig)

        return fig

    def _highlight_snp(
        self,
        ax: Any,
        fig: Any,
        snp_idx: int,
        n_snps: int,
        color: str,
    ) -> None:
        """Add visual highlight for a SNP's row/column in the heatmap.

        Draws rectangle borders around the row and column cells for the
        given SNP in the lower triangle.

        Args:
            ax: Axes object from backend.
            fig: Figure object from backend.
            snp_idx: Index of the SNP to highlight.
            n_snps: Total number of SNPs in the matrix.
            color: Highlight color.
        """
        # Compute all cell positions to highlight (x, y pairs)
        # Row cells: columns 0 to snp_idx, row = snp_idx
        row_cells = [(j, snp_idx) for j in range(snp_idx + 1)]
        # Column cells: column = snp_idx, rows snp_idx+1 to end (skip diagonal)
        col_cells = [(snp_idx, i) for i in range(snp_idx + 1, n_snps)]
        all_cells = row_cells + col_cells

        if self.backend_name == "matplotlib":
            from matplotlib.patches import Rectangle

            for x, y in all_cells:
                rect = Rectangle(
                    (x - 0.5, y - 0.5),
                    1.0,
                    1.0,
                    fill=False,
                    edgecolor=color,
                    linewidth=2,
                    zorder=10,
                )
                ax.add_patch(rect)

        elif self.backend_name == "plotly":
            for x, y in all_cells:
                fig.add_shape(
                    type="rect",
                    x0=x - 0.5,
                    x1=x + 0.5,
                    y0=y - 0.5,
                    y1=y + 0.5,
                    line=dict(color=color, width=2),
                    fillcolor="rgba(0,0,0,0)",
                )

        elif self.backend_name == "bokeh":
            for x, y in all_cells:
                ax.rect(
                    x=x,
                    y=y,
                    width=1,
                    height=1,
                    fill_alpha=0,
                    line_color=color,
                    line_width=2,
                )
