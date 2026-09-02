"""Miami figure: the request, the two mirrored panels, and the plan builder."""

from dataclasses import dataclass
from typing import Any, Optional, Tuple

from ._figure import FigurePlan, RegionHighlight
from ._manhattan_panel import ManhattanPanelSpec, manhattan_spec
from .backends.base import PlotBackend
from .backends.hover import HoverConfig
from .manhattan import PreparedManhattan


@dataclass(frozen=True)
class MiamiRequest:
    """One mirrored Manhattan figure, resolved by the plotter.

    ``top`` and ``bottom`` are prepared against one shared ``GenomeLayout``.
    ``rs_col`` names the id column the annotations index into and is set
    whenever either annotation tuple is non-empty.
    """

    top: PreparedManhattan
    bottom: PreparedManhattan
    hover: Optional[HoverConfig]
    rs_col: Optional[str]
    top_threshold: Optional[float]
    bottom_threshold: Optional[float]
    top_label: Optional[str]
    bottom_label: Optional[str]
    top_annotations: Tuple[str, ...]
    bottom_annotations: Tuple[str, ...]
    highlights: Tuple[Tuple[str, int, int], ...]
    highlight_color: str
    highlight_alpha: float
    figsize: Tuple[float, float]
    title: Optional[str]


@dataclass(frozen=True)
class MiamiPanel:
    """One half of a Miami figure: a Manhattan panel and its SNP annotations.

    ``rs_col`` is a column of the spec's frame whenever it is not None, and
    ``annotations`` are the ids in it to label.
    """

    spec: ManhattanPanelSpec
    rs_col: Optional[str]
    annotations: Tuple[str, ...]

    def draw(self, backend: PlotBackend, ax: Any) -> None:
        """Draw the Manhattan panel, then label the annotated SNPs."""
        self.spec.draw(backend, ax)
        if self.rs_col is None or not self.annotations:
            return
        frame = self.spec.prepared_df
        for _, row in frame[frame[self.rs_col].isin(self.annotations)].iterrows():
            backend.add_text(
                ax,
                x=row["_cumulative_pos"],
                y=row["neglog10p"],
                text=str(row[self.rs_col]),
                fontsize=8,
                ha="center",
                va="bottom",
            )


def miami_plan(req: MiamiRequest) -> FigurePlan:
    """Lay out the two mirrored panels and the highlights spanning both."""
    top = MiamiPanel(
        spec=manhattan_spec(
            req.top,
            significance_threshold=req.top_threshold,
            panel_label=req.top_label,
            hover=req.hover,
        ),
        rs_col=req.rs_col,
        annotations=req.top_annotations,
    )
    bottom = MiamiPanel(
        spec=manhattan_spec(
            req.bottom,
            significance_threshold=req.bottom_threshold,
            x_label="Chromosome",
            panel_label=req.bottom_label,
            panel_label_y_frac=0.05,
            invert_y=True,
            hover=req.hover,
        ),
        rs_col=req.rs_col,
        annotations=req.bottom_annotations,
    )
    offsets = req.top.layout.offsets
    highlights = [
        RegionHighlight(
            offsets[str(chrom)] + start,
            offsets[str(chrom)] + end,
            req.highlight_color,
            req.highlight_alpha,
        )
        for chrom, start, end in req.highlights
        if str(chrom) in offsets
    ]
    return FigurePlan(
        panels=[top, bottom],
        figsize=req.figsize,
        highlights=highlights,
        suptitle=req.title,
        top=0.92 if req.title else 0.95,
        hspace=0.05,
    )
