"""Tests for HoverDataBuilder hover data construction."""

import pandas as pd
import pytest

from pylocuszoom.backends.hover import (
    HoverConfig,
    HoverDataBuilder,
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
        """plotly_hovertemplate generates valid template string."""
        template = plotly_hovertemplate(hover_df)

        assert isinstance(template, str)
        assert template.endswith("<extra></extra>")

    def test_plotly_template_snp_first_bold(self, hover_df):
        """plotly_hovertemplate puts the first column in bold, unformatted."""
        template = plotly_hovertemplate(hover_df)

        assert template.startswith("<b>%{customdata[0]}</b>")

    def test_plotly_template_full_layout(self, hover_df):
        """Every column gets its name, index, and name-implied format."""
        assert plotly_hovertemplate(hover_df) == (
            "<b>%{customdata[0]}</b><br>"
            "Position: %{customdata[1]:,.0f}<br>"
            "P-value: %{customdata[2]:.2e}<br>"
            "R²: %{customdata[3]:.3f}<br>"
            "<extra></extra>"
        )

    def test_plotly_template_extra_cols_default_format(self):
        """A column with no recognised name renders without a format spec."""
        hover_df = pd.DataFrame({"SNP": ["rs123"], "Effect": [0.5]})

        assert plotly_hovertemplate(hover_df) == (
            "<b>%{customdata[0]}</b><br>Effect: %{customdata[1]}<br><extra></extra>"
        )

    def test_plotly_template_single_column(self):
        """A SNP-only frame still terminates the template correctly."""
        hover_df = pd.DataFrame({"SNP": ["rs123"]})

        assert plotly_hovertemplate(hover_df) == (
            "<b>%{customdata[0]}</b><br><extra></extra>"
        )


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

    def test_bokeh_tooltips_full_layout(self, hover_df):
        """Each column becomes a (display name, field reference) pair in order."""
        assert bokeh_tooltips(hover_df) == [
            ("SNP", "@{SNP}"),
            ("Position", "@{Position}{0,0}"),
            ("P-value", "@{P-value}{0.2e}"),
            ("R²", "@{R²}{0.3f}"),
        ]

    def test_bokeh_tooltips_key_prefix_namespaces_field_references(self, hover_df):
        """The prefix reaches the field reference but not the display name.

        BokehBackend.scatter relies on this to keep a hover column named "size"
        or "x" from shadowing the keys it sets for geometry and styling.
        """
        assert bokeh_tooltips(hover_df, key_prefix="hover_") == [
            ("SNP", "@{hover_SNP}"),
            ("Position", "@{hover_Position}{0,0}"),
            ("P-value", "@{hover_P-value}{0.2e}"),
            ("R²", "@{hover_R²}{0.3f}"),
        ]

    def test_bokeh_tooltips_extra_cols_no_format(self):
        """A column with no recognised name gets a bare field reference."""
        hover_df = pd.DataFrame({"SNP": ["rs123"], "Effect": [0.5]})

        assert bokeh_tooltips(hover_df) == [
            ("SNP", "@{SNP}"),
            ("Effect", "@{Effect}"),
        ]


class TestFormatDetection:
    """The column-name heuristic both backends share."""

    @pytest.mark.parametrize(
        "col_name,plotly_spec,bokeh_spec",
        [
            ("P-value", ".2e", "0.2e"),
            ("pval", ".2e", "0.2e"),
            ("p_value", ".2e", "0.2e"),
            ("PVAL", ".2e", "0.2e"),
            ("R²", ".3f", "0.3f"),
            ("r2", ".3f", "0.3f"),
            ("LD", ".3f", "0.3f"),
            ("ld_r2", ".3f", "0.3f"),
            ("Position", ",.0f", "0,0"),
            ("POS", ",.0f", "0,0"),
        ],
    )
    def test_recognised_names_get_their_format(self, col_name, plotly_spec, bokeh_spec):
        hover_df = pd.DataFrame({"SNP": ["rs123"], col_name: [1.0]})

        assert f":{plotly_spec}}}" in plotly_hovertemplate(hover_df)
        assert bokeh_tooltips(hover_df)[1] == (
            col_name,
            f"@{{{col_name}}}{{{bokeh_spec}}}",
        )

    @pytest.mark.parametrize(
        "col_name", ["Effect", "Gene", "beta", "MAF", "Consequence"]
    )
    def test_unrecognised_names_get_no_format(self, col_name):
        hover_df = pd.DataFrame({"SNP": ["rs123"], col_name: [1.0]})

        assert f"{col_name}: %{{customdata[1]}}" in plotly_hovertemplate(hover_df)
        assert bokeh_tooltips(hover_df)[1] == (col_name, f"@{{{col_name}}}")

    def test_p_value_name_match_is_exact_not_substring(self):
        """ "pval" formats as a p-value; "n_pval_tests" does not."""
        exact = pd.DataFrame({"SNP": ["rs1"], "pval": [1.0]})
        substring = pd.DataFrame({"SNP": ["rs1"], "n_pval_tests": [1.0]})

        assert ":.2e}" in plotly_hovertemplate(exact)
        assert ":.2e}" not in plotly_hovertemplate(substring)
