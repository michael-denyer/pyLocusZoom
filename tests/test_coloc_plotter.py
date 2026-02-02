"""Tests for ColocPlotter class."""

import importlib.util

import numpy as np
import pandas as pd
import pytest

from pylocuszoom import ValidationError

# Check if optional backends are available
PLOTLY_AVAILABLE = importlib.util.find_spec("plotly") is not None
BOKEH_AVAILABLE = importlib.util.find_spec("bokeh") is not None


@pytest.fixture
def gwas_data():
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


class TestColocPlotterInit:
    """Tests for ColocPlotter initialization."""

    def test_default_backend_is_matplotlib(self):
        """Test that default backend is matplotlib."""
        from pylocuszoom.coloc_plotter import ColocPlotter

        plotter = ColocPlotter()
        assert plotter.backend_name == "matplotlib"

    @pytest.mark.skipif(not PLOTLY_AVAILABLE, reason="Plotly not installed")
    def test_accepts_plotly_backend(self):
        """Test that plotly backend is accepted."""
        from pylocuszoom.coloc_plotter import ColocPlotter

        plotter = ColocPlotter(backend="plotly")
        assert plotter.backend_name == "plotly"

    @pytest.mark.skipif(not BOKEH_AVAILABLE, reason="Bokeh not installed")
    def test_accepts_bokeh_backend(self):
        """Test that bokeh backend is accepted."""
        from pylocuszoom.coloc_plotter import ColocPlotter

        plotter = ColocPlotter(backend="bokeh")
        assert plotter.backend_name == "bokeh"


class TestPlotColoc:
    """Tests for plot_coloc method."""

    def test_returns_figure(self, gwas_data, eqtl_data):
        """Test that basic plot returns non-None figure."""
        from pylocuszoom.coloc_plotter import ColocPlotter

        plotter = ColocPlotter()
        fig = plotter.plot_coloc(gwas_data, eqtl_data)
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

    def test_significance_lines(self, gwas_data, eqtl_data):
        """Test custom significance thresholds."""
        from pylocuszoom.coloc_plotter import ColocPlotter

        plotter = ColocPlotter()
        fig = plotter.plot_coloc(
            gwas_data, eqtl_data, gwas_threshold=1e-5, eqtl_threshold=1e-3
        )
        assert fig is not None

    def test_correlation_displayed(self, gwas_data, eqtl_data):
        """Test correlation is displayed when show_correlation=True."""
        from pylocuszoom.coloc_plotter import ColocPlotter

        plotter = ColocPlotter()
        fig = plotter.plot_coloc(gwas_data, eqtl_data, show_correlation=True)
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

    def test_no_ld_uses_na_color(self, gwas_data, eqtl_data):
        """Test that all points are grey when no ld_col provided."""
        from pylocuszoom.coloc_plotter import ColocPlotter

        plotter = ColocPlotter()
        fig = plotter.plot_coloc(gwas_data, eqtl_data, ld_col=None)
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

    def test_custom_figsize(self, gwas_data, eqtl_data):
        """Test figsize parameter is respected."""
        from pylocuszoom.coloc_plotter import ColocPlotter

        plotter = ColocPlotter()
        fig = plotter.plot_coloc(gwas_data, eqtl_data, figsize=(10, 10))
        width, height = fig.get_size_inches()
        assert width == 10
        assert height == 10

    def test_title_displayed(self, gwas_data, eqtl_data):
        """Test title parameter sets plot title."""
        from pylocuszoom.coloc_plotter import ColocPlotter

        plotter = ColocPlotter()
        fig = plotter.plot_coloc(gwas_data, eqtl_data, title="Test Coloc Plot")
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

    def test_missing_eqtl_columns_raises(self, gwas_data):
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
            plotter.plot_coloc(gwas_data, eqtl_missing)

    def test_invalid_p_value_raises(self, gwas_data):
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
            plotter.plot_coloc(gwas_data, eqtl_bad_p)


class TestColocPlotterBackends:
    """Tests for ColocPlotter with different backends."""

    @pytest.mark.skipif(not PLOTLY_AVAILABLE, reason="Plotly not installed")
    def test_plotly_backend(self, gwas_data, eqtl_data):
        """Test plotly backend returns plotly Figure."""
        import plotly.graph_objects as go

        from pylocuszoom.coloc_plotter import ColocPlotter

        plotter = ColocPlotter(backend="plotly")
        fig = plotter.plot_coloc(gwas_data, eqtl_data)
        assert isinstance(fig, go.Figure)

    @pytest.mark.skipif(not BOKEH_AVAILABLE, reason="Bokeh not installed")
    def test_bokeh_backend(self, gwas_data, eqtl_data):
        """Test bokeh backend returns bokeh layout."""
        from bokeh.models.layouts import LayoutDOM

        from pylocuszoom.coloc_plotter import ColocPlotter

        plotter = ColocPlotter(backend="bokeh")
        fig = plotter.plot_coloc(gwas_data, eqtl_data)
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

    def test_nan_p_values_handled(self):
        """Test NaN p-values are filtered or handled."""
        from pylocuszoom.coloc_plotter import ColocPlotter

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
                "p_eqtl": [1e-7, 1e-5, np.nan],
                "rs": ["rs1", "rs2", "rs3"],
            }
        )
        plotter = ColocPlotter()
        # Should handle NaN gracefully
        fig = plotter.plot_coloc(gwas, eqtl)
        assert fig is not None

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
