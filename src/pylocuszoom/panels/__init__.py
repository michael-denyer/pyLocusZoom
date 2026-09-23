"""One module per panel type: the value, the constructor, and its ``draw``.

Every panel resolves its mode, its hover contract, and its layout when it is
built, so each ``draw`` method issues backend primitives without probing the
frame again.  A plotter puts the panels it built on one
:class:`~.._figure.FigurePlan` and ``render_figure`` draws them in order.

The five regional panels are re-exported here because one plotter builds them
together and :data:`RegionalPanel` is the union over them.  The panels of the
other families have one caller each and are imported from their own module.
"""

from typing import List, Optional, Union

import pandas as pd

from ..config import LDHeatmapInput, PanelInputs, RegionConfig
from .association import AssociationInput, AssociationPanel
from .eqtl import EqtlPanel
from .finemapping import FinemappingPanel
from .genes import GenePanel
from .heatmap import HeatmapPanel

RegionalPanel = Union[
    AssociationPanel,
    FinemappingPanel,
    EqtlPanel,
    GenePanel,
    HeatmapPanel,
]


def optional_panels(
    inputs: PanelInputs,
    region: RegionConfig,
    *,
    genes_df: Optional[pd.DataFrame],
    exons_df: Optional[pd.DataFrame],
) -> List[RegionalPanel]:
    """Build the fine-mapping, eQTL and gene panels a figure asks for.

    They need no association panel, so a regional plot builds them, and
    validates their frames, before it runs PLINK for LD.

    Args:
        inputs: The caller's optional panel inputs.
        region: The figure's region.
        genes_df: Genes to draw, the caller's or fetched, or None for no track.
        exons_df: Exons for ``genes_df``, or None.

    Returns:
        The panels in figure order.
    """
    panels: List[RegionalPanel] = []
    if inputs.finemapping is not None:
        finemapping = inputs.finemapping
        panels.append(
            FinemappingPanel.from_frame(
                finemapping.data,
                region,
                finemapping.cs_col,
                chrom_col=finemapping.chrom_col,
            )
        )
    if inputs.eqtl is not None:
        eqtl = inputs.eqtl
        panels.append(
            EqtlPanel.from_frame(
                eqtl.data, region, eqtl.gene, eqtl.threshold, chrom_col=eqtl.chrom_col
            )
        )
    if genes_df is not None:
        panels.append(GenePanel.from_genes(genes_df, region, exons_df))
    return panels


def ld_heatmap_panels(
    heatmap: Optional[LDHeatmapInput],
    *,
    source: AssociationPanel,
    region: RegionConfig,
    association_height: float,
) -> List[HeatmapPanel]:
    """Build the LD heatmap panel under ``source``, if the figure asks for one.

    Args:
        heatmap: The caller's LD heatmap input, or None for no panel.
        source: The association panel whose SNPs place the heatmap.
        region: The figure's region.
        association_height: Height-ratio units of an association panel,
            which the heatmap's own height scales.

    Returns:
        No panel, or the one heatmap panel.
    """
    if heatmap is None:
        return []
    return [
        HeatmapPanel.from_matrix(
            heatmap.matrix,
            heatmap.snp_ids,
            source=source,
            region=region,
            height=association_height * heatmap.height,
            metric=heatmap.metric,
        )
    ]


__all__ = [
    "AssociationInput",
    "AssociationPanel",
    "EqtlPanel",
    "FinemappingPanel",
    "GenePanel",
    "HeatmapPanel",
    "RegionalPanel",
    "ld_heatmap_panels",
    "optional_panels",
]
