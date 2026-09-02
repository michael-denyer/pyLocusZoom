"""Tests for ColocPlotter class."""

import numpy as np
import pandas as pd
import pytest

from pylocuszoom import ValidationError


@pytest.fixture
def coloc_gwas_df():
    """Create sample GWAS data with 5 SNPs."""
    return pd.DataFrame(
        {
            "pos": [1000, 2000, 3000, 4000, 5000],
            "p_gwas": [1e-8, 1e-6, 1e-4, 0.01, 0.1],
            "rs": ["rs1", "rs2", "rs3", "rs4", "rs5"],
        }
    )


@pytest.fixture
def eqtl_data():
    """Create sample eQTL data with 5 overlapping SNPs."""
    return pd.DataFrame(
        {
            "pos": [1000, 2000, 3000, 4000, 5000],
            "p_eqtl": [1e-7, 1e-5, 1e-3, 0.02, 0.15],
            "rs": ["rs1", "rs2", "rs3", "rs4", "rs5"],
        }
    )


@pytest.fixture
def gwas_data_with_ld():
    """Create sample GWAS data with LD column."""
    return pd.DataFrame(
        {
            "pos": [1000, 2000, 3000, 4000, 5000],
            "p_gwas": [1e-8, 1e-6, 1e-4, 0.01, 0.1],
            "rs": ["rs1", "rs2", "rs3", "rs4", "rs5"],
            "ld": [1.0, 0.85, 0.65, 0.35, 0.15],
        }
    )


def _dashed_x(ax) -> list:
    """Return the x position of every vertical dashed line on a matplotlib axis."""
    return [
        line.get_xdata()[0]
        for line in ax.get_lines()
        if line.get_linestyle() == "--" and len(set(line.get_xdata())) == 1
    ]


def _dashed_y(ax) -> list:
    """Return the y position of every horizontal dashed line on a matplotlib axis."""
    return [
        line.get_ydata()[0]
        for line in ax.get_lines()
        if line.get_linestyle() == "--" and len(set(line.get_ydata())) == 1
    ]


class TestPlotColoc:
    """Tests for plot_coloc method."""

    def test_returns_figure(self, coloc_gwas_df, eqtl_data):
        """Test that basic plot returns non-None figure."""
        from pylocuszoom.coloc_plotter import ColocPlotter

        plotter = ColocPlotter()
        fig = plotter.plot_coloc(coloc_gwas_df, eqtl_data)
        assert fig is not None

    def test_ld_coloring(self, gwas_data_with_ld, eqtl_data):
        """Test that LD coloring works when ld_col is provided."""
        from pylocuszoom.coloc_plotter import ColocPlotter

        plotter = ColocPlotter()
        fig = plotter.plot_coloc(
            gwas_data_with_ld, eqtl_data, ld_col="ld", lead_snp="rs1"
        )
        assert fig is not None

    def test_lead_snp_highlighting(self, gwas_data_with_ld, eqtl_data):
        """Test that lead SNP is highlighted when specified."""
        from pylocuszoom.coloc_plotter import ColocPlotter

        plotter = ColocPlotter()
        fig = plotter.plot_coloc(
            gwas_data_with_ld, eqtl_data, ld_col="ld", lead_snp="rs1"
        )
        assert fig is not None

    def test_auto_select_lead_snp(self, gwas_data_with_ld, eqtl_data):
        """Test lead SNP auto-selection when ld_col provided but no lead_snp."""
        from pylocuszoom.coloc_plotter import ColocPlotter

        plotter = ColocPlotter()
        # Should auto-select rs1 (highest combined -log10p)
        fig = plotter.plot_coloc(gwas_data_with_ld, eqtl_data, ld_col="ld")
        assert fig is not None

    def test_significance_lines(self, coloc_gwas_df, eqtl_data):
        """Test custom significance thresholds."""
        from pylocuszoom.coloc_plotter import ColocPlotter

        plotter = ColocPlotter()
        fig = plotter.plot_coloc(
            coloc_gwas_df, eqtl_data, gwas_threshold=1e-5, eqtl_threshold=1e-3
        )
        assert fig is not None

    def test_thresholds_default_to_the_plotters_own(self, coloc_gwas_df, eqtl_data):
        """An unsupplied threshold inherits the plotter's, not a module default."""
        from pylocuszoom.coloc_plotter import ColocPlotter

        plotter = ColocPlotter(genomewide_threshold=1e-3, eqtl_threshold=1e-2)

        fig = plotter.plot_coloc(coloc_gwas_df, eqtl_data)

        ax = fig.get_axes()[0]
        assert _dashed_x(ax) == pytest.approx([3.0])
        assert _dashed_y(ax) == pytest.approx([2.0])

    def test_none_threshold_draws_no_line(self, coloc_gwas_df, eqtl_data):
        """None means no line, and does not fall back to the plotter default."""
        from pylocuszoom.coloc_plotter import ColocPlotter

        plotter = ColocPlotter()

        fig = plotter.plot_coloc(
            coloc_gwas_df, eqtl_data, gwas_threshold=None, eqtl_threshold=None
        )

        ax = fig.get_axes()[0]
        assert _dashed_x(ax) == []
        assert _dashed_y(ax) == []

    def test_correlation_displayed(self, coloc_gwas_df, eqtl_data):
        """Test correlation is displayed when show_correlation=True."""
        from pylocuszoom.coloc_plotter import ColocPlotter

        plotter = ColocPlotter()
        fig = plotter.plot_coloc(coloc_gwas_df, eqtl_data, show_correlation=True)
        assert fig is not None

    def test_correlation_skipped_small_n(self):
        """Test correlation is skipped when n < 3."""
        from pylocuszoom.coloc_plotter import ColocPlotter

        gwas_small = pd.DataFrame(
            {
                "pos": [1000, 2000],
                "p_gwas": [1e-8, 1e-6],
                "rs": ["rs1", "rs2"],
            }
        )
        eqtl_small = pd.DataFrame(
            {
                "pos": [1000, 2000],
                "p_eqtl": [1e-7, 1e-5],
                "rs": ["rs1", "rs2"],
            }
        )
        plotter = ColocPlotter()
        fig = plotter.plot_coloc(gwas_small, eqtl_small, show_correlation=True)
        assert fig is not None

    def test_no_ld_uses_na_color(self, coloc_gwas_df, eqtl_data):
        """Test that all points are grey when no ld_col provided."""
        from pylocuszoom.coloc_plotter import ColocPlotter

        plotter = ColocPlotter()
        fig = plotter.plot_coloc(coloc_gwas_df, eqtl_data, ld_col=None)
        assert fig is not None

    def test_custom_column_names(self):
        """Test non-default column names work."""
        from pylocuszoom.coloc_plotter import ColocPlotter

        gwas_custom = pd.DataFrame(
            {
                "position": [1000, 2000, 3000],
                "pval": [1e-8, 1e-6, 1e-4],
                "snp_id": ["rs1", "rs2", "rs3"],
            }
        )
        eqtl_custom = pd.DataFrame(
            {
                "position": [1000, 2000, 3000],
                "pval_eqtl": [1e-7, 1e-5, 1e-3],
                "snp_id": ["rs1", "rs2", "rs3"],
            }
        )
        plotter = ColocPlotter()
        fig = plotter.plot_coloc(
            gwas_custom,
            eqtl_custom,
            pos_col="position",
            gwas_p_col="pval",
            eqtl_p_col="pval_eqtl",
            rs_col="snp_id",
        )
        assert fig is not None

    def test_custom_figsize(self, coloc_gwas_df, eqtl_data):
        """Test figsize parameter is respected."""
        from pylocuszoom.coloc_plotter import ColocPlotter

        plotter = ColocPlotter()
        fig = plotter.plot_coloc(coloc_gwas_df, eqtl_data, figsize=(10, 10))
        width, height = fig.get_size_inches()
        assert width == 10
        assert height == 10

    def test_title_displayed(self, coloc_gwas_df, eqtl_data):
        """Test title parameter sets plot title."""
        from pylocuszoom.coloc_plotter import ColocPlotter

        plotter = ColocPlotter()
        fig = plotter.plot_coloc(coloc_gwas_df, eqtl_data, title="Test Coloc Plot")
        axes = fig.get_axes()
        assert len(axes) >= 1
        title_text = axes[0].get_title()
        assert "Test Coloc Plot" in title_text


class TestColocPlotterValidation:
    """Tests for ColocPlotter validation."""

    def test_empty_merge_raises_error(self):
        """Test that no overlapping positions raises ValueError."""
        from pylocuszoom.coloc_plotter import ColocPlotter

        gwas = pd.DataFrame(
            {
                "pos": [1000, 2000, 3000],
                "p_gwas": [1e-8, 1e-6, 1e-4],
                "rs": ["rs1", "rs2", "rs3"],
            }
        )
        eqtl = pd.DataFrame(
            {
                "pos": [4000, 5000, 6000],  # No overlap
                "p_eqtl": [1e-7, 1e-5, 1e-3],
                "rs": ["rs4", "rs5", "rs6"],
            }
        )
        plotter = ColocPlotter()
        with pytest.raises(ValueError, match="overlapping"):
            plotter.plot_coloc(gwas, eqtl)

    def test_missing_gwas_columns_raises(self, eqtl_data):
        """Test missing required columns in GWAS df raises ValidationError."""
        from pylocuszoom.coloc_plotter import ColocPlotter

        gwas_missing = pd.DataFrame(
            {
                "pos": [1000, 2000, 3000],
                # Missing p_gwas column
                "rs": ["rs1", "rs2", "rs3"],
            }
        )
        plotter = ColocPlotter()
        with pytest.raises(ValidationError, match="Missing columns"):
            plotter.plot_coloc(gwas_missing, eqtl_data)

    def test_missing_eqtl_columns_raises(self, coloc_gwas_df):
        """Test missing required columns in eQTL df raises ValidationError."""
        from pylocuszoom.coloc_plotter import ColocPlotter

        eqtl_missing = pd.DataFrame(
            {
                # Missing pos column
                "p_eqtl": [1e-7, 1e-5, 1e-3],
                "rs": ["rs1", "rs2", "rs3"],
            }
        )
        plotter = ColocPlotter()
        with pytest.raises(ValidationError, match="Missing columns"):
            plotter.plot_coloc(coloc_gwas_df, eqtl_missing)

    def test_invalid_p_value_raises(self, coloc_gwas_df):
        """Test p-values outside (0, 1] range raises ValidationError."""
        from pylocuszoom.coloc_plotter import ColocPlotter

        eqtl_bad_p = pd.DataFrame(
            {
                "pos": [1000, 2000, 3000],
                "p_eqtl": [0, 1e-5, 1e-3],  # 0 is invalid
                "rs": ["rs1", "rs2", "rs3"],
            }
        )
        plotter = ColocPlotter()
        with pytest.raises(ValidationError, match="values <= 0"):
            plotter.plot_coloc(coloc_gwas_df, eqtl_bad_p)


class TestColocPlotterBackends:
    """Tests for ColocPlotter with different backends."""

    def test_plotly_backend(self, coloc_gwas_df, eqtl_data):
        """Test plotly backend returns plotly Figure."""
        import plotly.graph_objects as go

        from pylocuszoom.coloc_plotter import ColocPlotter

        plotter = ColocPlotter(backend="plotly")
        fig = plotter.plot_coloc(coloc_gwas_df, eqtl_data)
        assert isinstance(fig, go.Figure)

    def test_bokeh_backend(self, coloc_gwas_df, eqtl_data):
        """Test bokeh backend returns bokeh layout."""
        from bokeh.models.layouts import LayoutDOM

        from pylocuszoom.coloc_plotter import ColocPlotter

        plotter = ColocPlotter(backend="bokeh")
        fig = plotter.plot_coloc(coloc_gwas_df, eqtl_data)
        assert isinstance(fig, LayoutDOM)

    def test_matplotlib_ld_legend(self, gwas_data_with_ld, eqtl_data):
        """Test matplotlib shows LD legend when ld_col provided."""
        from pylocuszoom.coloc_plotter import ColocPlotter

        plotter = ColocPlotter(backend="matplotlib")
        fig = plotter.plot_coloc(
            gwas_data_with_ld, eqtl_data, ld_col="ld", lead_snp="rs1"
        )
        # Legend should be present
        axes = fig.get_axes()
        assert len(axes) >= 1


class TestColocPlotterEdgeCases:
    """Tests for edge cases in ColocPlotter."""

    def test_single_overlapping_point(self):
        """Test works with single overlapping point (no correlation)."""
        from pylocuszoom.coloc_plotter import ColocPlotter

        gwas = pd.DataFrame(
            {
                "pos": [1000],
                "p_gwas": [1e-8],
                "rs": ["rs1"],
            }
        )
        eqtl = pd.DataFrame(
            {
                "pos": [1000],
                "p_eqtl": [1e-7],
                "rs": ["rs1"],
            }
        )
        plotter = ColocPlotter()
        fig = plotter.plot_coloc(gwas, eqtl)
        assert fig is not None

    def test_nan_p_values_rejected(self):
        """NaN p-values must raise ValidationError rather than blanking axes.

        See pyLocusZoom-3oy / pyLocusZoom-y0l.
        """
        from pylocuszoom.coloc_plotter import ColocPlotter
        from pylocuszoom.exceptions import ValidationError

        gwas = pd.DataFrame(
            {
                "pos": [1000, 2000, 3000],
                "p_gwas": [1e-8, np.nan, 1e-4],
                "rs": ["rs1", "rs2", "rs3"],
            }
        )
        eqtl = pd.DataFrame(
            {
                "pos": [1000, 2000, 3000],
                "p_eqtl": [1e-7, 1e-5, 1e-3],
                "rs": ["rs1", "rs2", "rs3"],
            }
        )
        plotter = ColocPlotter()
        with pytest.raises(ValidationError, match="null"):
            plotter.plot_coloc(gwas, eqtl)

    def test_lead_snp_not_in_data_raises(self, gwas_data_with_ld, eqtl_data):
        """Test invalid lead_snp raises ValueError."""
        from pylocuszoom.coloc_plotter import ColocPlotter

        plotter = ColocPlotter()
        with pytest.raises(ValueError, match="lead_snp"):
            plotter.plot_coloc(
                gwas_data_with_ld,
                eqtl_data,
                ld_col="ld",
                lead_snp="invalid_rs",
            )


@pytest.fixture
def gwas_data_with_effects():
    """GWAS data with effect sizes."""
    return pd.DataFrame(
        {
            "pos": [1000, 2000, 3000, 4000, 5000],
            "p_gwas": [1e-8, 1e-6, 0.01, 0.05, 0.1],
            "rs": ["rs1", "rs2", "rs3", "rs4", "rs5"],
            "beta_gwas": [
                0.5,
                -0.3,
                0.2,
                -0.1,
                0.4,
            ],  # positive, negative, positive, negative, positive
        }
    )


@pytest.fixture
def eqtl_data_with_effects():
    """eQTL data with effect sizes."""
    return pd.DataFrame(
        {
            "pos": [1000, 2000, 3000, 4000, 5000],
            "p_eqtl": [1e-7, 1e-5, 0.02, 0.04, 0.15],
            "rs": ["rs1", "rs2", "rs3", "rs4", "rs5"],
            "beta_eqtl": [
                0.4,
                -0.2,
                -0.1,
                -0.2,
                0.3,
            ],  # same, same, opposite, same, same
        }
    )


class TestEffectDirectionColoring:
    """Tests for effect direction coloring feature."""

    def test_color_by_effect_basic(
        self, gwas_data_with_effects, eqtl_data_with_effects
    ):
        """Test color_by_effect renders plot without error."""
        from pylocuszoom.coloc_plotter import ColocPlotter

        plotter = ColocPlotter()
        fig = plotter.plot_coloc(
            gwas_data_with_effects,
            eqtl_data_with_effects,
            color_by_effect=True,
            gwas_effect_col="beta_gwas",
            eqtl_effect_col="beta_eqtl",
        )
        assert fig is not None

    def test_effect_congruent_color(self):
        """Test same direction effects get green color."""
        from pylocuszoom.coloc_plotter import _get_effect_agreement_color
        from pylocuszoom.colors import EFFECT_CONGRUENT_COLOR

        # Both positive
        assert _get_effect_agreement_color(0.5, 0.4) == EFFECT_CONGRUENT_COLOR
        # Both negative
        assert _get_effect_agreement_color(-0.3, -0.2) == EFFECT_CONGRUENT_COLOR

    def test_effect_incongruent_color(self):
        """Test opposite direction effects get red color."""
        from pylocuszoom.coloc_plotter import _get_effect_agreement_color
        from pylocuszoom.colors import EFFECT_INCONGRUENT_COLOR

        # Positive GWAS, negative eQTL
        assert _get_effect_agreement_color(0.5, -0.2) == EFFECT_INCONGRUENT_COLOR
        # Negative GWAS, positive eQTL
        assert _get_effect_agreement_color(-0.3, 0.4) == EFFECT_INCONGRUENT_COLOR

    def test_color_by_effect_without_cols_raises(self, coloc_gwas_df, eqtl_data):
        """Test color_by_effect=True without effect cols raises ValueError."""
        from pylocuszoom.coloc_plotter import ColocPlotter

        plotter = ColocPlotter()

        # Missing both columns
        with pytest.raises(ValueError, match="color_by_effect.*requires"):
            plotter.plot_coloc(
                coloc_gwas_df,
                eqtl_data,
                color_by_effect=True,
            )

        # Missing eqtl_effect_col
        with pytest.raises(ValueError, match="color_by_effect.*requires"):
            plotter.plot_coloc(
                coloc_gwas_df,
                eqtl_data,
                color_by_effect=True,
                gwas_effect_col="beta_gwas",
            )

    def test_effect_nan_handled(self):
        """Test NaN effects get grey color."""
        from pylocuszoom.coloc_plotter import _get_effect_agreement_color
        from pylocuszoom.colors import LD_NA_COLOR

        # NaN GWAS effect
        assert _get_effect_agreement_color(np.nan, 0.4) == LD_NA_COLOR
        # NaN eQTL effect
        assert _get_effect_agreement_color(0.5, np.nan) == LD_NA_COLOR
        # Both NaN
        assert _get_effect_agreement_color(np.nan, np.nan) == LD_NA_COLOR

    def test_effect_legend_matplotlib(
        self, gwas_data_with_effects, eqtl_data_with_effects
    ):
        """Test effect direction legend shown for matplotlib backend."""
        from pylocuszoom.coloc_plotter import ColocPlotter

        plotter = ColocPlotter(backend="matplotlib")
        fig = plotter.plot_coloc(
            gwas_data_with_effects,
            eqtl_data_with_effects,
            color_by_effect=True,
            gwas_effect_col="beta_gwas",
            eqtl_effect_col="beta_eqtl",
        )
        # Check that legend exists on axis
        axes = fig.get_axes()
        assert len(axes) >= 1
        legend = axes[0].get_legend()
        assert legend is not None
        # Check legend has expected labels
        labels = [t.get_text() for t in legend.get_texts()]
        assert "Same direction" in labels
        assert "Opposite direction" in labels

    def test_missing_effect_column_raises(
        self, gwas_data_with_effects, eqtl_data_with_effects
    ):
        """Test missing effect column in data raises ValueError."""
        from pylocuszoom.coloc_plotter import ColocPlotter

        plotter = ColocPlotter()

        # Wrong GWAS effect column name
        with pytest.raises(ValueError, match="gwas_effect_col.*not found"):
            plotter.plot_coloc(
                gwas_data_with_effects,
                eqtl_data_with_effects,
                color_by_effect=True,
                gwas_effect_col="wrong_col",
                eqtl_effect_col="beta_eqtl",
            )


class TestH4PosteriorDisplay:
    """Tests for H4 posterior probability display."""

    def test_h4_displayed(self, coloc_gwas_df, eqtl_data):
        """Test H4 posterior is displayed when provided."""
        from pylocuszoom.coloc_plotter import ColocPlotter

        plotter = ColocPlotter()
        fig = plotter.plot_coloc(coloc_gwas_df, eqtl_data, h4_posterior=0.95)
        assert fig is not None

    def test_h4_range_validation(self, coloc_gwas_df, eqtl_data):
        """Test h4_posterior must be in [0, 1]."""
        from pylocuszoom.coloc_plotter import ColocPlotter

        plotter = ColocPlotter()

        # Invalid: < 0
        with pytest.raises(ValueError, match="h4_posterior"):
            plotter.plot_coloc(coloc_gwas_df, eqtl_data, h4_posterior=-0.1)

        # Invalid: > 1
        with pytest.raises(ValueError, match="h4_posterior"):
            plotter.plot_coloc(coloc_gwas_df, eqtl_data, h4_posterior=1.5)

    def test_h4_formatting(self, coloc_gwas_df, eqtl_data):
        """Test H4 is formatted to 3 decimal places."""
        from pylocuszoom.coloc_plotter import ColocPlotter

        plotter = ColocPlotter(backend="matplotlib")
        fig = plotter.plot_coloc(coloc_gwas_df, eqtl_data, h4_posterior=0.95123456)
        # Check that text annotation exists
        axes = fig.get_axes()
        assert len(axes) >= 1
        # Find H4 text in annotations
        texts = [child.get_text() for child in axes[0].texts]
        h4_texts = [t for t in texts if "H4 PP" in t]
        assert len(h4_texts) == 1
        assert "H4 PP = 0.951" in h4_texts[0]

    def test_h4_with_correlation(self, coloc_gwas_df, eqtl_data):
        """Test both H4 and correlation are displayed together."""
        from pylocuszoom.coloc_plotter import ColocPlotter

        plotter = ColocPlotter(backend="matplotlib")
        fig = plotter.plot_coloc(
            coloc_gwas_df, eqtl_data, h4_posterior=0.95, show_correlation=True
        )
        axes = fig.get_axes()
        texts = [child.get_text() for child in axes[0].texts]
        # Both annotations should be present
        h4_texts = [t for t in texts if "H4 PP" in t]
        corr_texts = [t for t in texts if "r = " in t]
        assert len(h4_texts) >= 1
        assert len(corr_texts) >= 1

    def test_h4_boundary_values(self, coloc_gwas_df, eqtl_data):
        """Test H4 boundary values are accepted."""
        from pylocuszoom.coloc_plotter import ColocPlotter

        plotter = ColocPlotter()

        # H4 = 0 (valid)
        fig = plotter.plot_coloc(coloc_gwas_df, eqtl_data, h4_posterior=0)
        assert fig is not None

        # H4 = 1 (valid)
        fig = plotter.plot_coloc(coloc_gwas_df, eqtl_data, h4_posterior=1)
        assert fig is not None


class TestColocConfigIntegration:
    """Tests for ColocConfig integration with ColocPlotter."""

    def test_config_to_plot_params(
        self, gwas_data_with_effects, eqtl_data_with_effects
    ):
        """Test ColocConfig values can drive plot_coloc parameters."""
        from pylocuszoom.coloc_plotter import ColocPlotter
        from pylocuszoom.config import ColocConfig

        config = ColocConfig(
            color_by_effect=True,
            gwas_effect_col="beta_gwas",
            eqtl_effect_col="beta_eqtl",
            h4_posterior=0.95,
            show_correlation=False,
            figsize=(10.0, 10.0),
        )

        plotter = ColocPlotter()
        fig = plotter.plot_coloc(
            gwas_data_with_effects,
            eqtl_data_with_effects,
            color_by_effect=config.color_by_effect,
            gwas_effect_col=config.gwas_effect_col,
            eqtl_effect_col=config.eqtl_effect_col,
            h4_posterior=config.h4_posterior,
            show_correlation=config.show_correlation,
            figsize=config.figsize,
        )
        assert fig is not None
        # Check figsize was applied
        width, height = fig.get_size_inches()
        assert width == 10.0
        assert height == 10.0

    def test_invalid_config_caught(self):
        """Test invalid ColocConfig raises ValidationError."""
        from pydantic import ValidationError

        from pylocuszoom.config import ColocConfig

        # color_by_effect without effect columns
        with pytest.raises(ValidationError, match="color_by_effect"):
            ColocConfig(color_by_effect=True)

        # h4_posterior out of range
        with pytest.raises(ValidationError, match="h4_posterior"):
            ColocConfig(h4_posterior=1.5)


class TestLeadSelectionRules:
    """When plot_coloc picks a lead variant, and when it picks none."""

    @staticmethod
    def _frames():
        gwas = pd.DataFrame(
            {
                "pos": [100, 200, 300],
                "p_gwas": [1e-8, 1e-3, 0.4],
                "rs": ["rs1", "rs2", "rs3"],
                "ld": [1.0, 0.6, 0.1],
            }
        )
        eqtl = pd.DataFrame(
            {
                "pos": [100, 200, 300],
                "p_eqtl": [1e-6, 1e-2, 0.5],
                "rs": ["rs1", "rs2", "rs3"],
            }
        )
        return gwas, eqtl

    def test_auto_selects_the_strongest_combined_signal_when_ld_is_present(self):
        from pylocuszoom.coloc_plotter import (
            ColocConfig,
            _merge_and_transform,
            _resolve_lead_idx,
        )

        gwas, eqtl = self._frames()
        config = ColocConfig(ld_col="ld")
        merged = _merge_and_transform(gwas, eqtl, config)

        lead = _resolve_lead_idx(merged, config)

        assert merged.data.loc[lead, merged.rs_col] == "rs1"

    def test_no_lead_without_ld_and_without_a_named_snp(self):
        from pylocuszoom.coloc_plotter import (
            ColocConfig,
            _merge_and_transform,
            _resolve_lead_idx,
        )

        gwas, eqtl = self._frames()
        config = ColocConfig(ld_col=None)
        merged = _merge_and_transform(gwas, eqtl, config)

        assert _resolve_lead_idx(merged, config) is None

    def test_a_named_lead_snp_beats_the_auto_selection(self):
        from pylocuszoom.coloc_plotter import (
            ColocConfig,
            _merge_and_transform,
            _resolve_lead_idx,
        )

        gwas, eqtl = self._frames()
        config = ColocConfig(ld_col="ld", lead_snp="rs3")
        merged = _merge_and_transform(gwas, eqtl, config)

        lead = _resolve_lead_idx(merged, config)

        assert merged.data.loc[lead, merged.rs_col] == "rs3"

    def test_named_lead_snp_without_an_rs_column_raises(self):
        from pylocuszoom.coloc_plotter import (
            ColocConfig,
            _merge_and_transform,
            _resolve_lead_idx,
        )

        gwas, eqtl = self._frames()
        gwas, eqtl = gwas.drop(columns=["rs"]), eqtl.drop(columns=["rs"])
        config = ColocConfig(rs_col=None, lead_snp="rs1")
        merged = _merge_and_transform(gwas, eqtl, config)

        with pytest.raises(ValueError, match="rs_col not found"):
            _resolve_lead_idx(merged, config)
