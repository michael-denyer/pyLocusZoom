"""Tests for LDHeatmapPlotter class."""

import numpy as np
import pandas as pd
import pytest

from pylocuszoom.ld_heatmap_plotter import LDHeatmapPlotter


@pytest.fixture
def sample_ld_matrix():
    """Create sample 5x5 LD matrix."""
    snp_ids = [f"rs{i}" for i in range(1, 6)]
    # Symmetric matrix with realistic LD pattern (high LD decays with distance)
    values = np.array(
        [
            [1.0, 0.9, 0.7, 0.4, 0.2],
            [0.9, 1.0, 0.8, 0.5, 0.3],
            [0.7, 0.8, 1.0, 0.6, 0.4],
            [0.4, 0.5, 0.6, 1.0, 0.7],
            [0.2, 0.3, 0.4, 0.7, 1.0],
        ]
    )
    return pd.DataFrame(values, index=snp_ids, columns=snp_ids)


@pytest.fixture
def matrix_with_nan():
    """Create LD matrix with missing values."""
    values = np.array(
        [
            [1.0, 0.8, np.nan],
            [0.8, 1.0, 0.5],
            [np.nan, 0.5, 1.0],
        ]
    )
    return pd.DataFrame(
        values, index=["rs1", "rs2", "rs3"], columns=["rs1", "rs2", "rs3"]
    )


@pytest.fixture
def small_ld_matrix():
    """Create small 3x3 LD matrix."""
    values = np.array(
        [
            [1.0, 0.8, 0.5],
            [0.8, 1.0, 0.6],
            [0.5, 0.6, 1.0],
        ]
    )
    return pd.DataFrame(
        values, index=["rs1", "rs2", "rs3"], columns=["rs1", "rs2", "rs3"]
    )


class TestLDHeatmapPlotterInit:
    """Tests for LDHeatmapPlotter initialization."""

    def test_default_backend_is_matplotlib(self):
        """Test that default backend is matplotlib."""
        plotter = LDHeatmapPlotter()
        assert plotter.backend_name == "matplotlib"

    def test_accepts_plotly_backend(self):
        """Test that plotly backend is accepted."""
        plotter = LDHeatmapPlotter(backend="plotly")
        assert plotter.backend_name == "plotly"

    def test_accepts_bokeh_backend(self):
        """Test that bokeh backend is accepted."""
        plotter = LDHeatmapPlotter(backend="bokeh")
        assert plotter.backend_name == "bokeh"

    def test_species_parameter(self):
        """Test that species parameter is stored."""
        plotter = LDHeatmapPlotter(species="human")
        assert plotter.species == "human"


class TestPlotLDHeatmap:
    """Tests for plot_ld_heatmap method."""

    def test_accepts_dataframe_matrix(self, sample_ld_matrix):
        """Test that DataFrame matrix is accepted."""
        plotter = LDHeatmapPlotter()
        fig = plotter.plot_ld_heatmap(sample_ld_matrix)
        assert fig is not None

    def test_accepts_numpy_array(self):
        """Test that numpy array is accepted."""
        plotter = LDHeatmapPlotter()
        matrix = np.array(
            [
                [1.0, 0.8, 0.5],
                [0.8, 1.0, 0.6],
                [0.5, 0.6, 1.0],
            ]
        )
        fig = plotter.plot_ld_heatmap(matrix)
        assert fig is not None

    def test_raises_for_non_square_matrix(self):
        """Test that non-square matrix raises ValueError."""
        plotter = LDHeatmapPlotter()
        matrix = np.array([[1.0, 0.5], [0.5, 1.0], [0.3, 0.4]])  # 3x2 matrix

        with pytest.raises(ValueError, match="square"):
            plotter.plot_ld_heatmap(matrix)

    def test_uses_dataframe_index_for_snp_ids(self, small_ld_matrix):
        """Test that DataFrame index is used for snp_ids when not provided."""
        plotter = LDHeatmapPlotter()
        fig = plotter.plot_ld_heatmap(small_ld_matrix)

        # Verify figure was created (index was used successfully)
        assert fig is not None
        axes = fig.get_axes()
        assert len(axes) >= 1

    def test_uses_provided_snp_ids(self):
        """Test that provided snp_ids are used."""
        plotter = LDHeatmapPlotter()
        matrix = np.array(
            [
                [1.0, 0.8],
                [0.8, 1.0],
            ]
        )
        custom_ids = ["snp_a", "snp_b"]
        fig = plotter.plot_ld_heatmap(matrix, snp_ids=custom_ids)
        assert fig is not None

    def test_lead_snp_parameter(self, small_ld_matrix):
        """Test that lead_snp parameter works."""
        plotter = LDHeatmapPlotter()
        fig = plotter.plot_ld_heatmap(small_ld_matrix, lead_snp="rs1")
        assert fig is not None

    def test_raises_for_invalid_lead_snp(self, small_ld_matrix):
        """Test that invalid lead_snp raises ValueError."""
        plotter = LDHeatmapPlotter()

        with pytest.raises(ValueError, match="lead_snp"):
            plotter.plot_ld_heatmap(small_ld_matrix, lead_snp="invalid_snp")

    def test_highlight_snps_parameter(self, sample_ld_matrix):
        """Test that highlight_snps parameter works."""
        plotter = LDHeatmapPlotter()
        fig = plotter.plot_ld_heatmap(
            sample_ld_matrix,
            lead_snp="rs1",
            highlight_snps=["rs2", "rs3"],
        )
        assert fig is not None

    def test_raises_for_invalid_highlight_snps(self, small_ld_matrix):
        """Test that invalid highlight_snps raises ValueError."""
        plotter = LDHeatmapPlotter()

        with pytest.raises(ValueError, match="highlight_snp"):
            plotter.plot_ld_heatmap(
                small_ld_matrix,
                highlight_snps=["invalid_snp"],
            )

    def test_metric_r2_label(self, small_ld_matrix):
        """Test that metric='r2' produces figure."""
        plotter = LDHeatmapPlotter()
        fig = plotter.plot_ld_heatmap(small_ld_matrix, metric="r2")
        assert fig is not None

    def test_metric_dprime_label(self, small_ld_matrix):
        """Test that metric='dprime' produces figure."""
        plotter = LDHeatmapPlotter()
        fig = plotter.plot_ld_heatmap(small_ld_matrix, metric="dprime")
        assert fig is not None

    def test_figsize_parameter(self, small_ld_matrix):
        """Test that figsize parameter is respected."""
        plotter = LDHeatmapPlotter()
        fig = plotter.plot_ld_heatmap(small_ld_matrix, figsize=(10, 10))

        width, height = fig.get_size_inches()
        assert width == 10
        assert height == 10

    def test_title_parameter(self, small_ld_matrix):
        """Test that title parameter is set."""
        plotter = LDHeatmapPlotter()
        fig = plotter.plot_ld_heatmap(small_ld_matrix, title="Test LD Heatmap")
        assert fig is not None

        # Verify title is set (matplotlib-specific)
        axes = fig.get_axes()
        assert len(axes) >= 1
        # Title is set on the axes
        title_text = axes[0].get_title()
        assert "Test LD Heatmap" in title_text

    def test_show_colorbar_true(self, small_ld_matrix):
        """Test that colorbar is shown when show_colorbar=True."""
        plotter = LDHeatmapPlotter()
        fig = plotter.plot_ld_heatmap(small_ld_matrix, show_colorbar=True)

        # Matplotlib figures with colorbar have more axes
        axes = fig.get_axes()
        assert len(axes) >= 2  # Main plot + colorbar

    def test_show_colorbar_false(self, small_ld_matrix):
        """Test that colorbar is hidden when show_colorbar=False."""
        plotter = LDHeatmapPlotter()
        fig = plotter.plot_ld_heatmap(small_ld_matrix, show_colorbar=False)

        # Without colorbar, should have only main plot axis
        axes = fig.get_axes()
        assert len(axes) == 1


class TestLDHeatmapBackends:
    """Tests for LDHeatmapPlotter with different backends."""

    def test_matplotlib_returns_figure(self, small_ld_matrix):
        """Test that matplotlib backend returns matplotlib Figure."""
        from matplotlib.figure import Figure

        plotter = LDHeatmapPlotter(backend="matplotlib")
        fig = plotter.plot_ld_heatmap(small_ld_matrix)

        assert isinstance(fig, Figure)

    def test_plotly_returns_figure(self, small_ld_matrix):
        """Test that plotly backend returns plotly Figure."""
        import plotly.graph_objects as go

        plotter = LDHeatmapPlotter(backend="plotly")
        fig = plotter.plot_ld_heatmap(small_ld_matrix)

        assert isinstance(fig, go.Figure)

    def test_bokeh_returns_figure(self, small_ld_matrix):
        """Test that bokeh backend returns bokeh layout."""
        from bokeh.models.layouts import LayoutDOM

        plotter = LDHeatmapPlotter(backend="bokeh")
        fig = plotter.plot_ld_heatmap(small_ld_matrix)

        assert isinstance(fig, LayoutDOM)

    def test_matplotlib_handles_nan_values(self, matrix_with_nan):
        """Test that matplotlib backend handles NaN values."""
        plotter = LDHeatmapPlotter(backend="matplotlib")
        fig = plotter.plot_ld_heatmap(matrix_with_nan)
        assert fig is not None

    def test_plotly_handles_nan_values(self, matrix_with_nan):
        """Test that plotly backend handles NaN values."""
        plotter = LDHeatmapPlotter(backend="plotly")
        fig = plotter.plot_ld_heatmap(matrix_with_nan)
        assert fig is not None

    def test_bokeh_handles_nan_values(self, matrix_with_nan):
        """Test that bokeh backend handles NaN values."""
        plotter = LDHeatmapPlotter(backend="bokeh")
        fig = plotter.plot_ld_heatmap(matrix_with_nan)
        assert fig is not None


class TestLDHeatmapEdgeCases:
    """Tests for edge cases in LDHeatmapPlotter."""

    def test_single_snp_matrix(self):
        """Test that single-SNP matrix is handled."""
        plotter = LDHeatmapPlotter()
        matrix = pd.DataFrame([[1.0]], index=["rs1"], columns=["rs1"])
        fig = plotter.plot_ld_heatmap(matrix)
        assert fig is not None

    def test_all_nan_matrix(self):
        """Test that all-NaN matrix is handled."""
        plotter = LDHeatmapPlotter()
        values = np.array(
            [
                [np.nan, np.nan],
                [np.nan, np.nan],
            ]
        )
        matrix = pd.DataFrame(values, index=["rs1", "rs2"], columns=["rs1", "rs2"])
        fig = plotter.plot_ld_heatmap(matrix)
        assert fig is not None

    def test_diagonal_only_values(self):
        """Test matrix with only diagonal values (identity matrix)."""
        plotter = LDHeatmapPlotter()
        matrix = pd.DataFrame(
            np.eye(3),
            index=["rs1", "rs2", "rs3"],
            columns=["rs1", "rs2", "rs3"],
        )
        fig = plotter.plot_ld_heatmap(matrix)
        assert fig is not None

    def test_empty_highlight_snps_list(self, small_ld_matrix):
        """Test that empty highlight_snps list works."""
        plotter = LDHeatmapPlotter()
        fig = plotter.plot_ld_heatmap(small_ld_matrix, highlight_snps=[])
        assert fig is not None

    def test_snp_ids_length_mismatch(self):
        """Test that mismatched snp_ids length raises ValueError."""
        plotter = LDHeatmapPlotter()
        matrix = np.array(
            [
                [1.0, 0.8, 0.5],
                [0.8, 1.0, 0.6],
                [0.5, 0.6, 1.0],
            ]
        )
        # Provide wrong number of snp_ids
        with pytest.raises(ValueError, match="length"):
            plotter.plot_ld_heatmap(
                matrix, snp_ids=["rs1", "rs2"]
            )  # Only 2 for 3x3 matrix


class TestLDHeatmapHighlighting:
    """Tests for SNP highlighting functionality."""

    def test_lead_snp_highlighting_matplotlib(self, sample_ld_matrix):
        """Test that lead SNP highlighting adds patches in matplotlib."""
        plotter = LDHeatmapPlotter(backend="matplotlib")
        fig = plotter.plot_ld_heatmap(sample_ld_matrix, lead_snp="rs3")

        axes = fig.get_axes()
        main_ax = axes[0]

        # Should have rectangle patches for highlighting
        patches = [p for p in main_ax.patches if hasattr(p, "get_edgecolor")]
        assert len(patches) > 0, "Lead SNP highlighting should add patches"

    def test_lead_snp_highlighting_plotly(self, sample_ld_matrix):
        """Test that lead SNP highlighting works in plotly."""
        plotter = LDHeatmapPlotter(backend="plotly")
        fig = plotter.plot_ld_heatmap(sample_ld_matrix, lead_snp="rs3")

        # Plotly should have shapes for highlighting
        assert fig is not None
        # Shapes are added for highlighting
        shapes = fig.layout.shapes
        assert shapes is not None and len(shapes) > 0

    def test_lead_snp_highlighting_bokeh(self, sample_ld_matrix):
        """Test that lead SNP highlighting works in bokeh."""
        plotter = LDHeatmapPlotter(backend="bokeh")
        fig = plotter.plot_ld_heatmap(sample_ld_matrix, lead_snp="rs3")

        # Just verify it completes without error
        assert fig is not None

    def test_multiple_highlights_matplotlib(self, sample_ld_matrix):
        """Test multiple SNP highlights with different colors."""
        plotter = LDHeatmapPlotter(backend="matplotlib")
        fig = plotter.plot_ld_heatmap(
            sample_ld_matrix,
            lead_snp="rs1",
            highlight_snps=["rs3", "rs5"],
        )

        axes = fig.get_axes()
        main_ax = axes[0]

        # Should have multiple patches for all highlights
        patches = [p for p in main_ax.patches if hasattr(p, "get_edgecolor")]
        assert len(patches) > 0
