"""Tests for HoverDataBuilder hover data construction."""

import pandas as pd
import pytest

from pylocuszoom.backends.hover import HoverConfig, HoverDataBuilder


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

    def test_build_dataframe_renames_columns(self, gwas_df, full_config):
        """build_dataframe returns DataFrame with standardized column names."""
        builder = HoverDataBuilder(full_config)
        hover_df = builder.build_dataframe(gwas_df)

        assert hover_df is not None
        assert "SNP" in hover_df.columns
        assert "Position" in hover_df.columns
        assert "P-value" in hover_df.columns
        assert "R²" in hover_df.columns
        # Original column names should not be present
        assert "rs" not in hover_df.columns
        assert "pvalue" not in hover_df.columns

    def test_build_dataframe_preserves_values(self, gwas_df, full_config):
        """build_dataframe preserves the actual data values."""
        builder = HoverDataBuilder(full_config)
        hover_df = builder.build_dataframe(gwas_df)

        assert list(hover_df["SNP"]) == ["rs123", "rs456", "rs789"]
        assert list(hover_df["Position"]) == [1000000, 2000000, 3000000]
        assert list(hover_df["P-value"]) == [1e-8, 1e-5, 0.05]
        assert list(hover_df["R²"]) == [1.0, 0.8, 0.2]

    def test_build_dataframe_skips_missing_columns(self, gwas_df):
        """build_dataframe skips columns not present in DataFrame."""
        config = HoverConfig(
            snp_col="nonexistent",  # Does not exist
            pos_col="position",
            p_col="pvalue",
        )
        builder = HoverDataBuilder(config)
        hover_df = builder.build_dataframe(gwas_df)

        assert hover_df is not None
        assert "SNP" not in hover_df.columns  # Skipped because column missing
        assert "Position" in hover_df.columns
        assert "P-value" in hover_df.columns

    def test_build_dataframe_returns_none_when_all_missing(self):
        """build_dataframe returns None when all configured columns are missing."""
        df = pd.DataFrame({"other": [1, 2, 3]})
        config = HoverConfig(
            snp_col="rs",
            pos_col="position",
            p_col="pvalue",
        )
        builder = HoverDataBuilder(config)
        hover_df = builder.build_dataframe(df)

        assert hover_df is None

    def test_build_dataframe_with_extra_cols(self, gwas_df):
        """build_dataframe includes extra columns with custom display names."""
        df = gwas_df.copy()
        df["beta"] = [0.5, -0.3, 0.1]
        df["maf"] = [0.25, 0.10, 0.45]

        config = HoverConfig(
            snp_col="rs",
            p_col="pvalue",
            extra_cols={"beta": "Effect", "maf": "MAF"},
        )
        builder = HoverDataBuilder(config)
        hover_df = builder.build_dataframe(df)

        assert "Effect" in hover_df.columns
        assert "MAF" in hover_df.columns
        assert list(hover_df["Effect"]) == [0.5, -0.3, 0.1]
        assert list(hover_df["MAF"]) == [0.25, 0.10, 0.45]

    def test_build_dataframe_maintains_column_order(self, gwas_df, full_config):
        """build_dataframe returns columns in consistent order: SNP, Position, P-value, R², extras."""
        builder = HoverDataBuilder(full_config)
        hover_df = builder.build_dataframe(gwas_df)

        columns = list(hover_df.columns)
        assert columns == ["SNP", "Position", "P-value", "R²"]

    def test_build_dataframe_partial_config(self, gwas_df):
        """build_dataframe works with partial configuration."""
        config = HoverConfig(snp_col="rs", p_col="pvalue")
        builder = HoverDataBuilder(config)
        hover_df = builder.build_dataframe(gwas_df)

        assert list(hover_df.columns) == ["SNP", "P-value"]


class TestPlotlyTemplateGeneration:
    """Tests for Plotly hovertemplate generation."""

    @pytest.fixture
    def hover_df(self):
        """Sample hover DataFrame."""
        return pd.DataFrame(
            {
                "SNP": ["rs123"],
                "Position": [1000000],
                "P-value": [1e-8],
                "R²": [0.85],
            }
        )

    def test_plotly_template_basic_structure(self, hover_df):
        """build_plotly_template generates valid template string."""
        config = HoverConfig(
            snp_col="rs", pos_col="position", p_col="pvalue", ld_col="r2"
        )
        builder = HoverDataBuilder(config)
        template = builder.build_plotly_template(hover_df)

        assert isinstance(template, str)
        assert "<extra></extra>" in template

    def test_plotly_template_snp_first_bold(self, hover_df):
        """build_plotly_template puts SNP first in bold."""
        config = HoverConfig(snp_col="rs", pos_col="position")
        builder = HoverDataBuilder(config)
        template = builder.build_plotly_template(hover_df)

        assert template.startswith("<b>%{customdata[0]}</b>")

    def test_plotly_template_pvalue_scientific_notation(self, hover_df):
        """build_plotly_template uses .2e format for P-value."""
        config = HoverConfig(snp_col="rs", p_col="pvalue")
        builder = HoverDataBuilder(config)
        template = builder.build_plotly_template(hover_df)

        # P-value should use .2e format
        assert ":.2e}" in template

    def test_plotly_template_position_comma_format(self, hover_df):
        """build_plotly_template uses comma format for Position."""
        config = HoverConfig(snp_col="rs", pos_col="position")
        builder = HoverDataBuilder(config)
        template = builder.build_plotly_template(hover_df)

        # Position should use comma format
        assert ":,.0f}" in template

    def test_plotly_template_r2_decimal_format(self, hover_df):
        """build_plotly_template uses .3f format for R-squared."""
        config = HoverConfig(snp_col="rs", ld_col="r2")
        builder = HoverDataBuilder(config)
        template = builder.build_plotly_template(hover_df)

        # R² should use .3f format
        assert ":.3f}" in template

    def test_plotly_template_customdata_indices(self, hover_df):
        """build_plotly_template uses correct customdata indices."""
        config = HoverConfig(
            snp_col="rs", pos_col="position", p_col="pvalue", ld_col="r2"
        )
        builder = HoverDataBuilder(config)
        template = builder.build_plotly_template(hover_df)

        # Each column should get sequential index
        assert "customdata[0]" in template  # SNP
        assert "customdata[1]" in template  # Position
        assert "customdata[2]" in template  # P-value
        assert "customdata[3]" in template  # R²

    def test_plotly_template_extra_cols_default_format(self):
        """build_plotly_template uses default format for extra columns."""
        hover_df = pd.DataFrame({"SNP": ["rs123"], "Effect": [0.5]})
        config = HoverConfig(snp_col="rs", extra_cols={"beta": "Effect"})
        builder = HoverDataBuilder(config)
        template = builder.build_plotly_template(hover_df)

        # Extra column should NOT have format specifier (just value)
        # Should be customdata[1] without format
        assert "%{customdata[1]}" in template
        # Should NOT have .2e or .3f
        assert template.count(":.") == 0 or "customdata[1]}" in template


class TestBokehTooltipsGeneration:
    """Tests for Bokeh tooltips list generation."""

    @pytest.fixture
    def hover_df(self):
        """Sample hover DataFrame."""
        return pd.DataFrame(
            {
                "SNP": ["rs123"],
                "Position": [1000000],
                "P-value": [1e-8],
                "R²": [0.85],
            }
        )

    def test_bokeh_tooltips_returns_list_of_tuples(self, hover_df):
        """build_bokeh_tooltips returns list of (name, format) tuples."""
        config = HoverConfig(
            snp_col="rs", pos_col="position", p_col="pvalue", ld_col="r2"
        )
        builder = HoverDataBuilder(config)
        tooltips = builder.build_bokeh_tooltips(hover_df)

        assert isinstance(tooltips, list)
        assert all(isinstance(t, tuple) and len(t) == 2 for t in tooltips)

    def test_bokeh_tooltips_pvalue_scientific_notation(self, hover_df):
        """build_bokeh_tooltips uses {0.2e} for P-value."""
        config = HoverConfig(snp_col="rs", p_col="pvalue")
        builder = HoverDataBuilder(config)
        tooltips = builder.build_bokeh_tooltips(hover_df)

        pvalue_tooltip = next((t for t in tooltips if t[0] == "P-value"), None)
        assert pvalue_tooltip is not None
        assert "{0.2e}" in pvalue_tooltip[1]

    def test_bokeh_tooltips_position_comma_format(self, hover_df):
        """build_bokeh_tooltips uses {0,0} for Position."""
        config = HoverConfig(snp_col="rs", pos_col="position")
        builder = HoverDataBuilder(config)
        tooltips = builder.build_bokeh_tooltips(hover_df)

        pos_tooltip = next((t for t in tooltips if t[0] == "Position"), None)
        assert pos_tooltip is not None
        assert "{0,0}" in pos_tooltip[1]

    def test_bokeh_tooltips_r2_decimal_format(self, hover_df):
        """build_bokeh_tooltips uses {0.3f} for R-squared."""
        config = HoverConfig(snp_col="rs", ld_col="r2")
        builder = HoverDataBuilder(config)
        tooltips = builder.build_bokeh_tooltips(hover_df)

        r2_tooltip = next((t for t in tooltips if "R" in t[0] or "r" in t[0]), None)
        assert r2_tooltip is not None
        assert "{0.3f}" in r2_tooltip[1]

    def test_bokeh_tooltips_column_reference(self, hover_df):
        """build_bokeh_tooltips uses @{column} format for references."""
        config = HoverConfig(snp_col="rs", p_col="pvalue")
        builder = HoverDataBuilder(config)
        tooltips = builder.build_bokeh_tooltips(hover_df)

        # SNP should reference @{SNP}
        snp_tooltip = next((t for t in tooltips if t[0] == "SNP"), None)
        assert snp_tooltip is not None
        assert "@{SNP}" in snp_tooltip[1] or "@SNP" in snp_tooltip[1]

    def test_bokeh_tooltips_extra_cols_no_format(self):
        """build_bokeh_tooltips uses no format for extra columns."""
        hover_df = pd.DataFrame({"SNP": ["rs123"], "Effect": [0.5]})
        config = HoverConfig(snp_col="rs", extra_cols={"beta": "Effect"})
        builder = HoverDataBuilder(config)
        tooltips = builder.build_bokeh_tooltips(hover_df)

        effect_tooltip = next((t for t in tooltips if t[0] == "Effect"), None)
        assert effect_tooltip is not None
        # Should be @Effect or @{Effect} without format specifier
        assert "@" in effect_tooltip[1]
        assert "{0." not in effect_tooltip[1]

    def test_bokeh_tooltips_preserves_order(self, hover_df):
        """build_bokeh_tooltips preserves column order from hover_df."""
        config = HoverConfig(
            snp_col="rs", pos_col="position", p_col="pvalue", ld_col="r2"
        )
        builder = HoverDataBuilder(config)
        tooltips = builder.build_bokeh_tooltips(hover_df)

        names = [t[0] for t in tooltips]
        assert names == ["SNP", "Position", "P-value", "R²"]


class TestFormatSpecs:
    """Tests for format specification constants."""

    def test_format_specs_exist(self):
        """HoverDataBuilder has FORMAT_SPECS class attribute."""
        assert hasattr(HoverDataBuilder, "FORMAT_SPECS")

    def test_format_specs_pvalue(self):
        """FORMAT_SPECS has p_value key with .2e format."""
        assert "p_value" in HoverDataBuilder.FORMAT_SPECS
        assert HoverDataBuilder.FORMAT_SPECS["p_value"] == ".2e"

    def test_format_specs_r2(self):
        """FORMAT_SPECS has r2 key with .3f format."""
        assert "r2" in HoverDataBuilder.FORMAT_SPECS
        assert HoverDataBuilder.FORMAT_SPECS["r2"] == ".3f"

    def test_format_specs_position(self):
        """FORMAT_SPECS has position key with ,.0f format."""
        assert "position" in HoverDataBuilder.FORMAT_SPECS
        assert HoverDataBuilder.FORMAT_SPECS["position"] == ",.0f"
