"""Colocalization scatter plot for GWAS-eQTL visualization.

Creates scatter plots comparing GWAS -log10(p) vs eQTL -log10(p)
with points colored by LD to the lead SNP.
"""

from typing import Any, Optional

from ._figure import FigurePlan, render_figure
from ._plotter_utils import (
    DEFAULT_EQTL_THRESHOLD,
    DEFAULT_GENOMEWIDE_THRESHOLD,
    UNSET,
    ThresholdArg,
    resolve_threshold,
)
from .backends import BackendType, get_backend
from .config import ColocConfig
from .exceptions import ValidationError
from .panels.coloc import ColocPanel
from .utils import DataFrameLike, to_pandas


class ColocPlotter:
    """Colocalization scatter plot generator.

    Creates scatter plots comparing GWAS -log10(p) vs eQTL -log10(p)
    with points colored by LD to the lead SNP.

    Supports multiple rendering backends:
    - matplotlib (default): Static publication-quality plots
    - plotly: Interactive HTML with hover tooltips
    - bokeh: Interactive HTML for dashboards

    Args:
        backend: Plotting backend ('matplotlib', 'plotly', or 'bokeh').
        genomewide_threshold: P-value threshold for the GWAS significance line.
        eqtl_threshold: P-value threshold for the eQTL significance line.

    Example:
        >>> plotter = ColocPlotter()
        >>> fig = plotter.plot_coloc(
        ...     gwas_df, eqtl_df, config=ColocConfig(lead_snp="rs12345")
        ... )
        >>> fig.savefig("coloc.png", dpi=150)
    """

    def __init__(
        self,
        backend: BackendType = "matplotlib",
        genomewide_threshold: float = DEFAULT_GENOMEWIDE_THRESHOLD,
        eqtl_threshold: float = DEFAULT_EQTL_THRESHOLD,
    ):
        """Initialize the colocalization plotter."""
        self._backend = get_backend(backend)
        self.genomewide_threshold = genomewide_threshold
        self.eqtl_threshold = eqtl_threshold

    def plot_coloc(
        self,
        gwas_df: DataFrameLike,
        eqtl_df: DataFrameLike,
        *,
        config: ColocConfig = ColocConfig(),
        gwas_threshold: ThresholdArg = UNSET,
        eqtl_threshold: ThresholdArg = UNSET,
        title: Optional[str] = None,
    ) -> Any:
        """Create GWAS-eQTL colocalization scatter plot.

        Args:
            gwas_df: GWAS results DataFrame with positions and p-values.
            eqtl_df: eQTL results DataFrame with positions and p-values.
            config: :class:`~pylocuszoom.ColocConfig` naming the columns of
                both frames, the lead SNP, the colouring and annotations, and
                the figure size.
            gwas_threshold: Significance threshold for the GWAS line. Defaults
                to the plotter's ``genomewide_threshold``; pass None to draw no
                line.
            eqtl_threshold: Significance threshold for the eQTL line. Defaults
                to the plotter's ``eqtl_threshold``; pass None to draw no line.
            title: Plot title.

        Returns:
            Figure object (type depends on backend).

        Raises:
            ValidationError: If required or named columns are missing or
                invalid, no position is in both frames, ``config.lead_snp``
                is not found, or a threshold is outside (0, 1].

        Example:
            >>> from pylocuszoom import ColocConfig
            >>> fig = plotter.plot_coloc(
            ...     gwas_df,
            ...     eqtl_df,
            ...     config=ColocConfig(ld_col="ld", lead_snp="rs12345"),
            ... )
            >>> # With effect coloring
            >>> fig = plotter.plot_coloc(
            ...     gwas_df,
            ...     eqtl_df,
            ...     config=ColocConfig(
            ...         color_by_effect=True,
            ...         gwas_effect_col="beta_gwas",
            ...         eqtl_effect_col="beta_eqtl",
            ...     ),
            ... )
        """
        thresholds = {
            "gwas_threshold": resolve_threshold(
                gwas_threshold, self.genomewide_threshold
            ),
            "eqtl_threshold": resolve_threshold(eqtl_threshold, self.eqtl_threshold),
        }
        for name, value in thresholds.items():
            if value is not None and not 0 < value <= 1:
                raise ValidationError(f"{name} must be in (0, 1], got {value}")
        panel = ColocPanel.from_frames(
            to_pandas(gwas_df), to_pandas(eqtl_df), config, title=title, **thresholds
        )
        return render_figure(
            self._backend, FigurePlan(panels=[panel], figsize=config.figsize)
        )
