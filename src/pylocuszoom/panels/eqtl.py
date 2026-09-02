"""The regional eQTL panel: effect-signed markers over the same region."""

from dataclasses import dataclass
from typing import Any, Optional

import pandas as pd

from .._plotter_utils import add_significance_line
from ..backends.base import PlotBackend
from ..backends.composition import LegendEntry, eqtl_legend_entries
from ..backends.hover import HoverConfig, HoverDataBuilder
from ..colors import EQTL_MARKER_COLOR, get_eqtl_color
from ..config import RegionConfig
from ..eqtl import prepare_eqtl_for_plotting
from ._shared import REGIONAL_LINE_ALPHA


@dataclass(frozen=True)
class EqtlPanel:
    """Prepared eQTL panel; ``gene`` is the gene the data was filtered to.

    ``effect_col`` is None when the frame carries no effect sizes, which is
    the single-colour diamond mode.
    """

    data: pd.DataFrame
    height: float
    gene: Optional[str]
    threshold: float
    effect_col: Optional[str]
    hover: HoverConfig

    @classmethod
    def from_frame(
        cls,
        df: pd.DataFrame,
        region: RegionConfig,
        gene: Optional[str],
        threshold: float,
    ) -> "EqtlPanel":
        """Validate, gene- and region-filter, and transform raw eQTL results.

        Raises:
            EQTLValidationError: If required columns are missing, or ``gene``
                is given and the frame has no ``gene`` column.
        """
        data = prepare_eqtl_for_plotting(
            df, gene=gene, chrom=region.chrom, start=region.start, end=region.end
        )
        extra_cols = {
            col: label
            for col, label in (("effect_size", "Effect"), ("gene", "Gene"))
            if col in data.columns
        }
        return cls(
            data=data,
            height=2.0,
            gene=gene,
            threshold=threshold,
            effect_col="effect_size" if "effect_size" in data.columns else None,
            hover=HoverConfig(
                pos_col="pos" if "pos" in data.columns else None,
                p_col="p_value" if "p_value" in data.columns else None,
                extra_cols=extra_cols,
            ),
        )

    def draw(self, backend: PlotBackend, ax: Any) -> None:
        """Draw eQTL points, their legend, and the eQTL significance line."""
        data = self.data
        if not data.empty:
            hover_builder = HoverDataBuilder(self.hover)
            if self.effect_col is not None:
                effect = data[self.effect_col]
                for subset, marker in (
                    (data[effect >= 0], "^"),
                    (data[effect < 0], "v"),
                ):
                    if not subset.empty:
                        backend.scatter(
                            ax,
                            subset["pos"],
                            subset["neglog10p"],
                            colors=subset[self.effect_col]
                            .apply(get_eqtl_color)
                            .tolist(),
                            sizes=50,
                            marker=marker,
                            edgecolor="black",
                            linewidth=0.5,
                            zorder=2,
                            hover_data=hover_builder.build_dataframe(subset),
                        )
                backend.add_legend(
                    ax,
                    eqtl_legend_entries(),
                    loc="upper right",
                    title="eQTL effect",
                )
            else:
                label = f"eQTL ({self.gene})" if self.gene else "eQTL"
                backend.scatter(
                    ax,
                    data["pos"],
                    data["neglog10p"],
                    colors=EQTL_MARKER_COLOR,
                    sizes=60,
                    marker="D",
                    edgecolor="black",
                    linewidth=0.5,
                    zorder=2,
                    hover_data=hover_builder.build_dataframe(data),
                )
                backend.add_legend(
                    ax,
                    [LegendEntry(label, EQTL_MARKER_COLOR, marker="D")],
                    loc="upper right",
                )
        backend.set_ylabel(ax, r"$-\log_{10}$ P (eQTL)")
        add_significance_line(backend, ax, self.threshold, alpha=REGIONAL_LINE_ALPHA)
