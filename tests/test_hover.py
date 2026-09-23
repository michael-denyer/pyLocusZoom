"""Tests for HoverDataBuilder hover data construction."""

import pandas as pd
import pytest

from pylocuszoom.backends.hover import (
    HoverConfig,
    HoverData,
    HoverDataBuilder,
    HoverRole,
    bokeh_tooltips,
    plotly_hovertemplate,
)


class TestHoverConfig:
    """Tests for HoverConfig dataclass."""

    def test_default_values(self):
        """HoverConfig has correct defaults."""
        config = HoverConfig()
        assert config.snp_col is None
        assert config.pos_col is None
        assert config.p_col is None
        assert config.ld_col is None
        assert config.extra_cols == {}

    def test_with_all_columns(self):
        """HoverConfig accepts all column mappings."""
        config = HoverConfig(
            snp_col="rs",
            pos_col="position",
            p_col="pvalue",
            ld_col="r2",
            extra_cols={"beta": "Effect"},
        )
        assert config.snp_col == "rs"
        assert config.pos_col == "position"
        assert config.p_col == "pvalue"
        assert config.ld_col == "r2"
        assert config.extra_cols == {"beta": "Effect"}


class TestHoverDataBuilder:
    """Tests for HoverDataBuilder hover data construction."""

    @pytest.fixture
    def gwas_df(self):
        """Sample GWAS DataFrame for testing."""
        return pd.DataFrame(
            {
                "rs": ["rs123", "rs456", "rs789"],
                "position": [1000000, 2000000, 3000000],
                "pvalue": [1e-8, 1e-5, 0.05],
                "r2": [1.0, 0.8, 0.2],
            }
        )

    @pytest.fixture
    def full_config(self):
        """Config with all standard columns mapped."""
        return HoverConfig(
            snp_col="rs",
            pos_col="position",
            p_col="pvalue",
            ld_col="r2",
        )

    def test_build_renames_columns(self, gwas_df, full_config):
        """build returns DataFrame with standardized column names."""
        builder = HoverDataBuilder(full_config)
        hover_df = builder.build(gwas_df).frame

        assert hover_df is not None
        assert "SNP" in hover_df.columns
        assert "Position" in hover_df.columns
        assert "P-value" in hover_df.columns
        assert "R²" in hover_df.columns
        # Original column names should not be present
        assert "rs" not in hover_df.columns
        assert "pvalue" not in hover_df.columns

    def test_build_preserves_values(self, gwas_df, full_config):
        """build preserves the actual data values."""
        builder = HoverDataBuilder(full_config)
        hover_df = builder.build(gwas_df).frame

        assert list(hover_df["SNP"]) == ["rs123", "rs456", "rs789"]
        assert list(hover_df["Position"]) == [1000000, 2000000, 3000000]
        assert list(hover_df["P-value"]) == [1e-8, 1e-5, 0.05]
        assert list(hover_df["R²"]) == [1.0, 0.8, 0.2]

    def test_build_skips_missing_columns(self, gwas_df):
        """build skips columns not present in DataFrame."""
        config = HoverConfig(
            snp_col="nonexistent",  # Does not exist
            pos_col="position",
            p_col="pvalue",
        )
        builder = HoverDataBuilder(config)
        hover_df = builder.build(gwas_df).frame

        assert hover_df is not None
        assert "SNP" not in hover_df.columns  # Skipped because column missing
        assert "Position" in hover_df.columns
        assert "P-value" in hover_df.columns

    def test_build_returns_none_when_all_missing(self):
        """build returns None when all configured columns are missing."""
        df = pd.DataFrame({"other": [1, 2, 3]})
        config = HoverConfig(
            snp_col="rs",
            pos_col="position",
            p_col="pvalue",
        )
        builder = HoverDataBuilder(config)
        hover_df = builder.build(df)

        assert hover_df is None

    def test_build_with_extra_cols(self, gwas_df):
        """build includes extra columns with custom display names."""
        df = gwas_df.copy()
        df["beta"] = [0.5, -0.3, 0.1]
        df["maf"] = [0.25, 0.10, 0.45]

        config = HoverConfig(
            snp_col="rs",
            p_col="pvalue",
            extra_cols={"beta": "Effect", "maf": "MAF"},
        )
        builder = HoverDataBuilder(config)
        hover_df = builder.build(df).frame

        assert "Effect" in hover_df.columns
        assert "MAF" in hover_df.columns
        assert list(hover_df["Effect"]) == [0.5, -0.3, 0.1]
        assert list(hover_df["MAF"]) == [0.25, 0.10, 0.45]

    def test_build_maintains_column_order(self, gwas_df, full_config):
        """build returns columns in consistent order: SNP, Position, P-value, R², extras."""
        builder = HoverDataBuilder(full_config)
        hover_df = builder.build(gwas_df).frame

        columns = list(hover_df.columns)
        assert columns == ["SNP", "Position", "P-value", "R²"]

    def test_build_records_each_columns_role(self, gwas_df):
        """Each column's role comes from the config field it was mapped by."""
        df = gwas_df.assign(beta=[0.5, -0.3, 0.1])
        config = HoverConfig(
            snp_col="rs",
            pos_col="position",
            p_col="pvalue",
            ld_col="r2",
            extra_cols={"beta": "Effect"},
        )

        hover = HoverDataBuilder(config).build(df)

        assert hover.roles == (
            HoverRole.ID,
            HoverRole.POSITION,
            HoverRole.P_VALUE,
            HoverRole.R2,
            HoverRole.PLAIN,
        )

    def test_build_partial_config(self, gwas_df):
        """build works with partial configuration."""
        config = HoverConfig(snp_col="rs", p_col="pvalue")
        builder = HoverDataBuilder(config)
        hover_df = builder.build(gwas_df).frame

        assert list(hover_df.columns) == ["SNP", "P-value"]


FULL_HOVER = HoverData(
    pd.DataFrame(
        {"SNP": ["rs123"], "Position": [1000000], "P-value": [1e-8], "R²": [0.85]}
    ),
    (HoverRole.ID, HoverRole.POSITION, HoverRole.P_VALUE, HoverRole.R2),
)


class TestPlotlyTemplateGeneration:
    """Tests for Plotly hovertemplate generation."""

    def test_plotly_template_full_layout(self):
        """Every column gets its name, index and role format; the id is bold."""
        assert plotly_hovertemplate(FULL_HOVER) == (
            "<b>SNP: %{customdata[0]}</b><br>"
            "Position: %{customdata[1]:,.0f}<br>"
            "P-value: %{customdata[2]:.2e}<br>"
            "R²: %{customdata[3]:.3f}<br>"
            "<extra></extra>"
        )

    def test_the_first_column_is_not_assumed_to_be_the_id(self):
        """A frame with no SNP id labels and formats its first column."""
        hover = HoverData(
            pd.DataFrame({"Position": [1000000], "PIP": [0.9]}),
            (HoverRole.POSITION, HoverRole.PLAIN),
        )

        assert plotly_hovertemplate(hover) == (
            "Position: %{customdata[0]:,.0f}<br>PIP: %{customdata[1]}<br>"
            "<extra></extra>"
        )


class TestBokehTooltipsGeneration:
    """Tests for Bokeh tooltips list generation."""

    def test_bokeh_tooltips_full_layout(self):
        """Each column becomes a (display name, field reference) pair in order."""
        assert bokeh_tooltips(FULL_HOVER) == [
            ("SNP", "@{SNP}"),
            ("Position", "@{Position}{0,0}"),
            ("P-value", "@{P-value}{0.2e}"),
            ("R²", "@{R²}{0.3f}"),
        ]

    def test_bokeh_tooltips_key_prefix_namespaces_field_references(self):
        """The prefix reaches the field reference but not the display name.

        BokehBackend.scatter relies on this to keep a hover column named "size"
        or "x" from shadowing the keys it sets for geometry and styling.
        """
        assert bokeh_tooltips(FULL_HOVER, key_prefix="hover_") == [
            ("SNP", "@{hover_SNP}"),
            ("Position", "@{hover_Position}{0,0}"),
            ("P-value", "@{hover_P-value}{0.2e}"),
            ("R²", "@{hover_R²}{0.3f}"),
        ]


class TestFormatsFollowRoles:
    """A column's format depends on its role, never on its display name."""

    @pytest.mark.parametrize("name", ["Exposure", "Fold change", "pval", "LD"])
    def test_a_plain_column_is_unformatted_whatever_its_name(self, name):
        hover = HoverData(pd.DataFrame({name: [1.0]}), (HoverRole.PLAIN,))

        assert plotly_hovertemplate(hover) == (
            f"{name}: %{{customdata[0]}}<br><extra></extra>"
        )
        assert bokeh_tooltips(hover) == [(name, f"@{{{name}}}")]
