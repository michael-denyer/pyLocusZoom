"""One module per panel type: the value, the constructor, and its ``draw``.

Every panel resolves its mode, its hover contract, and its layout when it is
built, so each ``draw`` method issues backend primitives without probing the
frame again.  A plotter puts the panels it built on one
:class:`~.._figure.FigurePlan` and ``render_figure`` draws them in order.

The five regional panels are re-exported here because one plotter builds them
together and :data:`RegionalPanel` is the union over them.  The panels of the
other families have one caller each and are imported from their own module.
"""

from typing import Union

from .association import AssociationPanel, hover_for_association
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

__all__ = [
    "AssociationPanel",
    "EqtlPanel",
    "FinemappingPanel",
    "GenePanel",
    "HeatmapPanel",
    "RegionalPanel",
    "hover_for_association",
]
