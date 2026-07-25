"""Hover data and tooltip construction for the interactive backends.

``HoverDataBuilder`` turns a caller's column mapping into a display-named hover
DataFrame. ``plotly_hovertemplate`` and ``bokeh_tooltips`` turn that DataFrame
into each backend's native tooltip spec, so the column-name-to-number-format
heuristic is owned here rather than re-derived inside every ``scatter``.
"""

from dataclasses import dataclass, field
from typing import List, Optional, Tuple

import pandas as pd

# Display-name substrings mapped to the number format each backend wants.
# Checked in order: an exact p-value name, then any r^2/LD hint, then position.
_PLOTLY_FORMATS = {"p_value": ".2e", "r2": ".3f", "position": ",.0f"}
_BOKEH_FORMATS = {"p_value": "0.2e", "r2": "0.3f", "position": "0,0"}

_P_VALUE_NAMES = ("p-value", "pval", "p_value")
_R2_HINTS = ("r2", "r²", "ld")


def _format_key(col_name: str) -> Optional[str]:
    """Classify a hover column by display name, or None for default format."""
    col_lower = col_name.lower()
    if col_lower in _P_VALUE_NAMES:
        return "p_value"
    if any(hint in col_lower for hint in _R2_HINTS):
        return "r2"
    if "pos" in col_lower:
        return "position"
    return None


def plotly_hovertemplate(hover_df: pd.DataFrame) -> str:
    """Build a Plotly hovertemplate over ``hover_df`` as customdata.

    The first column is rendered bold and unformatted; it is the SNP identifier
    by ``HoverDataBuilder`` convention. Later columns get the format their
    display name implies.

    Args:
        hover_df: DataFrame returned by ``HoverDataBuilder.build_dataframe``.

    Returns:
        Plotly hovertemplate string referencing ``customdata`` by position.
    """
    parts = []
    for i, col in enumerate(hover_df.columns):
        if i == 0:
            parts.append(f"<b>%{{customdata[{i}]}}</b>")
            continue
        key = _format_key(col)
        spec = f":{_PLOTLY_FORMATS[key]}" if key else ""
        parts.append(f"{col}: %{{customdata[{i}]{spec}}}")
    parts.append("<extra></extra>")
    return "<br>".join(parts)


def bokeh_tooltips(
    hover_df: pd.DataFrame, key_prefix: str = ""
) -> List[Tuple[str, str]]:
    """Build Bokeh ``HoverTool`` tooltips over ``hover_df``.

    Args:
        hover_df: DataFrame returned by ``HoverDataBuilder.build_dataframe``.
        key_prefix: Prefix applied to each ``ColumnDataSource`` key, letting a
            caller namespace hover columns away from its own keys.

    Returns:
        List of ``(display_name, field_reference)`` tuples.
    """
    tooltips = []
    for col in hover_df.columns:
        key = _format_key(col)
        spec = f"{{{_BOKEH_FORMATS[key]}}}" if key else ""
        tooltips.append((col, f"@{{{key_prefix}{col}}}{spec}"))
    return tooltips


@dataclass
class HoverConfig:
    """Configuration for hover data column mapping.

    Maps source DataFrame column names to standardized display names for tooltips.

    Attributes:
        snp_col: Column name for SNP identifiers (displayed as "SNP").
        pos_col: Column name for genomic position (displayed as "Position").
        p_col: Column name for p-value (displayed as "P-value").
        ld_col: Column name for LD/R-squared (displayed as "R²").
        extra_cols: Additional columns to include, mapping source name to display name.
    """

    snp_col: Optional[str] = None
    pos_col: Optional[str] = None
    p_col: Optional[str] = None
    ld_col: Optional[str] = None
    extra_cols: dict[str, str] = field(default_factory=dict)


class HoverDataBuilder:
    """Builder for the display-named hover DataFrame a backend renders.

    Holds one ``HoverConfig`` so a caller can build hover data for several
    frames (all points, then the lead SNP) under the same column mapping.
    """

    # Standard column mappings (source config attr -> display name)
    _COLUMN_MAPPING = {
        "snp_col": "SNP",
        "pos_col": "Position",
        "p_col": "P-value",
        "ld_col": "R²",
    }

    def __init__(self, config: HoverConfig) -> None:
        """Initialize builder with column configuration.

        Args:
            config: HoverConfig with column name mappings.
        """
        self.config = config

    def build_dataframe(self, df: pd.DataFrame) -> Optional[pd.DataFrame]:
        """Build standardized hover DataFrame with renamed columns.

        Extracts configured columns from the input DataFrame, renames them to
        standardized display names, and returns a new DataFrame. Columns that
        don't exist in the input are skipped gracefully.

        Args:
            df: Input DataFrame containing hover data columns.

        Returns:
            DataFrame with renamed columns, or None if no configured columns exist.
        """
        result_data = {}

        # Process standard columns in order
        for config_attr, display_name in self._COLUMN_MAPPING.items():
            source_col = getattr(self.config, config_attr)
            if source_col is not None and source_col in df.columns:
                result_data[display_name] = df[source_col].values

        # Process extra columns
        for source_col, display_name in self.config.extra_cols.items():
            if source_col in df.columns:
                result_data[display_name] = df[source_col].values

        if not result_data:
            return None

        return pd.DataFrame(result_data)
