"""Hover data and tooltip construction for the interactive backends.

``HoverDataBuilder`` turns a caller's column mapping into ``HoverData``: the
display-named columns and the role each one plays. ``plotly_hovertemplate``
and ``bokeh_tooltips`` format each column by its role, so both backends show
the same fields in the same formats and neither guesses from a display name.
"""

from dataclasses import dataclass, field
from enum import Enum
from typing import List, Optional, Tuple

import pandas as pd


class HoverRole(Enum):
    """What a hover column holds, which decides how it is formatted."""

    ID = "id"
    POSITION = "position"
    P_VALUE = "p_value"
    R2 = "r2"
    PLAIN = "plain"


_PLOTLY_FORMATS = {
    HoverRole.POSITION: ":,.0f",
    HoverRole.P_VALUE: ":.2e",
    HoverRole.R2: ":.3f",
}
_BOKEH_FORMATS = {
    HoverRole.POSITION: "{0,0}",
    HoverRole.P_VALUE: "{0.2e}",
    HoverRole.R2: "{0.3f}",
}


@dataclass(frozen=True)
class HoverData:
    """The columns a tooltip shows, under their display names, with their roles.

    Attributes:
        frame: One column per tooltip field, in display order.
        roles: The role of each column of ``frame``, in the same order.
    """

    frame: pd.DataFrame
    roles: Tuple[HoverRole, ...]


def plotly_hovertemplate(hover: HoverData) -> str:
    """Build a Plotly hovertemplate over ``hover.frame`` as customdata.

    Every field is labelled with its display name; the SNP id line is bold.

    Args:
        hover: Columns and roles from ``HoverDataBuilder.build``.

    Returns:
        Plotly hovertemplate string referencing ``customdata`` by position.
    """
    parts = []
    for i, (col, role) in enumerate(zip(hover.frame.columns, hover.roles)):
        line = f"{col}: %{{customdata[{i}]{_PLOTLY_FORMATS.get(role, '')}}}"
        parts.append(f"<b>{line}</b>" if role is HoverRole.ID else line)
    parts.append("<extra></extra>")
    return "<br>".join(parts)


def bokeh_tooltips(hover: HoverData, key_prefix: str = "") -> List[Tuple[str, str]]:
    """Build Bokeh ``HoverTool`` tooltips over ``hover.frame``.

    Args:
        hover: Columns and roles from ``HoverDataBuilder.build``.
        key_prefix: Prefix applied to each ``ColumnDataSource`` key, letting a
            caller namespace hover columns away from its own keys.

    Returns:
        List of ``(display_name, field_reference)`` tuples.
    """
    return [
        (col, f"@{{{key_prefix}{col}}}{_BOKEH_FORMATS.get(role, '')}")
        for col, role in zip(hover.frame.columns, hover.roles)
    ]


@dataclass
class HoverConfig:
    """Configuration for hover data column mapping.

    Maps source DataFrame column names to standardized display names for tooltips.

    Attributes:
        snp_col: Column name for SNP identifiers (displayed as "SNP").
        pos_col: Column name for genomic position (displayed as "Position").
        p_col: Column name for p-value (displayed as "P-value").
        ld_col: Column name for LD/R-squared (displayed as "R²").
        extra_cols: Additional columns to include, mapping source name to
            display name. They are shown unformatted.
    """

    snp_col: Optional[str] = None
    pos_col: Optional[str] = None
    p_col: Optional[str] = None
    ld_col: Optional[str] = None
    extra_cols: dict[str, str] = field(default_factory=dict)


class HoverDataBuilder:
    """Builder for the hover columns and roles a backend renders.

    Holds one ``HoverConfig`` so a caller can build hover data for several
    frames (all points, then the lead SNP) under the same column mapping.
    """

    # Standard column mappings: config attr -> (display name, role)
    _COLUMN_MAPPING = {
        "snp_col": ("SNP", HoverRole.ID),
        "pos_col": ("Position", HoverRole.POSITION),
        "p_col": ("P-value", HoverRole.P_VALUE),
        "ld_col": ("R²", HoverRole.R2),
    }

    def __init__(self, config: HoverConfig) -> None:
        """Initialize builder with column configuration.

        Args:
            config: HoverConfig with column name mappings.
        """
        self.config = config

    def build(self, df: pd.DataFrame) -> Optional[HoverData]:
        """Build the display-named hover columns and their roles.

        Extracts configured columns from the input DataFrame, renames them to
        standardized display names, and records each one's role. Columns that
        don't exist in the input are skipped.

        Args:
            df: Input DataFrame containing hover data columns.

        Returns:
            The hover columns and roles, or None if no configured column exists.
        """
        columns = {}
        roles = []
        standard = (
            (getattr(self.config, attr), name, role)
            for attr, (name, role) in self._COLUMN_MAPPING.items()
        )
        extra = (
            (source, name, HoverRole.PLAIN)
            for source, name in self.config.extra_cols.items()
        )
        for source_col, display_name, role in (*standard, *extra):
            if source_col is not None and source_col in df.columns:
                columns[display_name] = df[source_col].values
                roles.append(role)

        if not columns:
            return None

        return HoverData(pd.DataFrame(columns), tuple(roles))
