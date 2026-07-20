"""Hover data builder for consistent hover tooltip generation.

Provides a unified interface for building hover data across Plotly and Bokeh backends.
"""

from dataclasses import dataclass, field
from typing import Optional

import pandas as pd


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
    """Builder for constructing hover data and templates across backends.

    Provides consistent hover tooltip generation for Plotly and Bokeh backends
    with automatic format detection for common column types.

    Attributes:
        FORMAT_SPECS: Format specifiers for different data types.
    """

    FORMAT_SPECS = {
        "p_value": ".2e",
        "r2": ".3f",
        "position": ",.0f",
    }

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

    def build_plotly_template(self, hover_df: pd.DataFrame) -> str:
        """Generate Plotly hovertemplate string.

        Creates a template string with appropriate format specifiers for each
        column type. SNP appears first in bold, followed by other columns with
        their values.

        Args:
            hover_df: DataFrame returned by build_dataframe().

        Returns:
            Plotly hovertemplate string with customdata references.
        """
        columns = hover_df.columns.tolist()
        parts = []

        for i, col in enumerate(columns):
            if i == 0:
                # First column (SNP) in bold
                parts.append(f"<b>%{{customdata[{i}]}}</b>")
            else:
                fmt = self._detect_plotly_format(col)
                if fmt:
                    parts.append(f"{col}: %{{customdata[{i}]:{fmt}}}")
                else:
                    parts.append(f"{col}: %{{customdata[{i}]}}")

        parts.append("<extra></extra>")
        return "<br>".join(parts)

    def build_bokeh_tooltips(self, hover_df: pd.DataFrame) -> list[tuple[str, str]]:
        """Generate Bokeh tooltips list.

        Creates a list of (name, format) tuples for Bokeh HoverTool configuration.
        Each tuple contains the display name and the Bokeh format string with
        appropriate specifiers.

        Args:
            hover_df: DataFrame returned by build_dataframe().

        Returns:
            List of (name, format_string) tuples for Bokeh tooltips.
        """
        tooltips = []

        for col in hover_df.columns:
            fmt = self._detect_bokeh_format(col)
            if fmt:
                tooltips.append((col, f"@{{{col}}}{{{fmt}}}"))
            else:
                tooltips.append((col, f"@{col}"))

        return tooltips

    def _detect_plotly_format(self, col_name: str) -> Optional[str]:
        """Detect appropriate Plotly format specifier for column.

        Uses heuristics based on column name to determine formatting:
        - P-value columns: scientific notation (.2e)
        - R²/LD columns: 3 decimal places (.3f)
        - Position columns: comma-separated (.0f)

        Args:
            col_name: Display name of the column.

        Returns:
            Plotly format specifier string, or None for default format.
        """
        col_lower = col_name.lower()

        if col_lower in ("p-value", "pval", "p_value"):
            return self.FORMAT_SPECS["p_value"]
        elif any(x in col_lower for x in ("r2", "r²", "ld")):
            return self.FORMAT_SPECS["r2"]
        elif "pos" in col_lower:
            return self.FORMAT_SPECS["position"]

        return None

    def _detect_bokeh_format(self, col_name: str) -> Optional[str]:
        """Detect appropriate Bokeh format code for column.

        Uses heuristics based on column name to determine formatting:
        - P-value columns: scientific notation (0.2e)
        - R²/LD columns: 3 decimal places (0.3f)
        - Position columns: comma-separated (0,0)

        Args:
            col_name: Display name of the column.

        Returns:
            Bokeh format code string, or None for default format.
        """
        col_lower = col_name.lower()

        if col_lower in ("p-value", "pval", "p_value"):
            return "0.2e"
        elif any(x in col_lower for x in ("r2", "r²", "ld")):
            return "0.3f"
        elif "pos" in col_lower:
            return "0,0"

        return None
