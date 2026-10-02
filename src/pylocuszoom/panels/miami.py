"""Miami figure: the request, the two mirrored panels, and the plan builder."""

import warnings
from dataclasses import dataclass
from typing import Any, Optional, Tuple, Union

from .._figure import FigurePlan, RegionHighlight
from ..backends.base import PlotBackend
from ..backends.hover import HoverConfig
from ..config import GenomeWideStyle
from ..manhattan import PreparedManhattan
from ..utils import normalize_chrom
from .manhattan import ManhattanPanelSpec, share_y_max


@dataclass(frozen=True)
class MiamiRequest:
    """One mirrored Manhattan figure, resolved by the plotter.

    ``top`` and ``bottom`` are prepared against one shared ``GenomeLayout``.
    ``rs_col`` names the id column the annotations index into and is set
    whenever either annotation tuple is non-empty. Each ``highlights`` entry
    is ``(chrom, start, end)`` with ``1 <= start <= end``.
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
    highlights: Tuple[Tuple[Union[int, str], int, int], ...]
    highlight_color: str
    highlight_alpha: float
    figsize: Tuple[float, float]
    title: Optional[str]
    style: GenomeWideStyle = GenomeWideStyle()


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
        frame, x_col = self.spec.prepared.frame, self.spec.prepared.x_col
        for _, row in frame[frame[self.rs_col].isin(self.annotations)].iterrows():
            backend.add_text(
                ax,
                x=row[x_col],
                y=row["neglog10p"],
                text=str(row[self.rs_col]),
                fontsize=8,
                ha="center",
                va="bottom",
            )


def miami_plan(req: MiamiRequest) -> FigurePlan:
    """Lay out the two mirrored panels and the highlights spanning both.

    A highlight stops at the last plotted position of its chromosome, where
    that chromosome's extent on the shared axis ends. One with nothing plotted
    under it (no data on the chromosome, or a start past its last plotted
    position) is skipped with a ``UserWarning`` naming the region.
    """
    top_spec, bottom_spec = share_y_max(
        [
            ManhattanPanelSpec(
                req.top,
                significance_threshold=req.top_threshold,
                panel_label=req.top_label,
                hover=req.hover,
                style=req.style,
            ),
            ManhattanPanelSpec(
                req.bottom,
                significance_threshold=req.bottom_threshold,
                x_label="Chromosome",
                panel_label=req.bottom_label,
                panel_label_y_frac=0.05,
                invert_y=True,
                hover=req.hover,
                style=req.style,
            ),
        ],
        req.style.y_headroom,
    )
    top = MiamiPanel(spec=top_spec, rs_col=req.rs_col, annotations=req.top_annotations)
    bottom = MiamiPanel(
        spec=bottom_spec, rs_col=req.rs_col, annotations=req.bottom_annotations
    )
    layout = req.top.layout
    highlights = []
    # A plain loop: the stacklevel below counts ``miami_plan`` and
    # ``plot_miami``, and a comprehension adds a frame on Python 3.10 and 3.11.
    for chrom, start, end in req.highlights:
        name = normalize_chrom(chrom)
        max_pos = layout.max_positions.get(name)
        if max_pos is not None and start <= max_pos:
            highlights.append(
                RegionHighlight(
                    layout.offsets[name] + start,
                    layout.offsets[name] + min(end, max_pos),
                    req.highlight_color,
                    req.highlight_alpha,
                )
            )
            continue
        if max_pos is None:
            why = f"chromosome {name} has no plotted data"
        else:
            why = (
                f"it starts past the last plotted position on chromosome {name} "
                f"({max_pos})"
            )
        warnings.warn(
            f"Highlight region {chrom}:{start}-{end} skipped; {why}",
            UserWarning,
            stacklevel=3,
        )
    return FigurePlan(
        panels=[top, bottom],
        figsize=req.figsize,
        highlights=highlights,
        suptitle=req.title,
        title_fontsize=req.style.title_fontsize,
        title_fontweight=req.style.title_fontweight,
        top=0.92 if req.title else 0.95,
        hspace=0.05,
    )
