"""Tests for LDHeatmapPlotter class."""

import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
import pytest

from pylocuszoom.backends import BUILTIN_BACKENDS
from pylocuszoom.colors import LEAD_SNP_HIGHLIGHT_COLOR, SECONDARY_HIGHLIGHT_COLOR
from pylocuszoom.ld_heatmap_plotter import LDHeatmapPlotter
from tests.conftest import FIGURE_TYPES

LEAD_OUTLINE = LEAD_SNP_HIGHLIGHT_COLOR.lower()
SECONDARY_OUTLINE = SECONDARY_HIGHLIGHT_COLOR.lower()


def drawn_cells(fig):
    """Return the heatmap values drawn on the main axis, NaN where masked."""
    (image,) = fig.get_axes()[0].images
    return np.ma.filled(image.get_array().astype(float), np.nan)


def outline_colors(fig):
    """Return the hex edge colour of every highlight rectangle, in draw order."""
    return [
        mcolors.to_hex(patch.get_edgecolor()) for patch in fig.get_axes()[0].patches
    ]


def tick_labels(fig):
    """Return the x and y tick label text of the main axis."""
    ax = fig.get_axes()[0]
    return (
        [tick.get_text() for tick in ax.get_xticklabels()],
        [tick.get_text() for tick in ax.get_yticklabels()],
    )


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


class TestPlotLDHeatmap:
    """Tests for plot_ld_heatmap method."""

    @pytest.mark.parametrize("as_dataframe", [True, False])
    def test_dataframe_and_array_draw_the_same_cells(
        self, small_ld_matrix, as_dataframe
    ):
        """Assert a DataFrame and its values draw the same lower triangle and labels."""
        plotter = LDHeatmapPlotter()
        matrix = small_ld_matrix if as_dataframe else small_ld_matrix.values

        fig = plotter.plot_ld_heatmap(matrix, snp_ids=["rs1", "rs2", "rs3"])

        np.testing.assert_array_equal(
            drawn_cells(fig),
            np.array(
                [
                    [1.0, np.nan, np.nan],
                    [0.8, 1.0, np.nan],
                    [0.5, 0.6, 1.0],
                ]
            ),
        )
        assert tick_labels(fig) == (["rs1", "rs2", "rs3"], ["rs1", "rs2", "rs3"])

    def test_raises_for_non_square_matrix(self):
        """Test that non-square matrix raises ValueError."""
        plotter = LDHeatmapPlotter()
        matrix = np.array([[1.0, 0.5], [0.5, 1.0], [0.3, 0.4]])  # 3x2 matrix

        with pytest.raises(ValueError, match="square"):
            plotter.plot_ld_heatmap(matrix)

    def test_uses_dataframe_index_for_snp_ids(self, small_ld_matrix):
        """Assert the DataFrame index labels the axes when snp_ids is not given."""
        plotter = LDHeatmapPlotter()
        fig = plotter.plot_ld_heatmap(small_ld_matrix)

        expected = list(small_ld_matrix.index)
        assert tick_labels(fig) == (expected, expected)

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

        assert tick_labels(fig) == (custom_ids, custom_ids)

    def test_lead_snp_parameter(self, small_ld_matrix):
        """Assert the lead SNP is outlined in the lead colour, one cell per SNP."""
        plotter = LDHeatmapPlotter()
        fig = plotter.plot_ld_heatmap(small_ld_matrix, lead_snp="rs1")

        assert outline_colors(fig) == [LEAD_OUTLINE] * 3

    def test_raises_for_invalid_lead_snp(self, small_ld_matrix):
        """Test that invalid lead_snp raises ValueError."""
        plotter = LDHeatmapPlotter()

        with pytest.raises(ValueError, match="lead_snp"):
            plotter.plot_ld_heatmap(small_ld_matrix, lead_snp="invalid_snp")

    def test_highlight_snps_parameter(self, sample_ld_matrix):
        """Assert every highlighted SNP is outlined in the secondary colour."""
        plotter = LDHeatmapPlotter()
        fig = plotter.plot_ld_heatmap(
            sample_ld_matrix,
            lead_snp="rs1",
            highlight_snps=["rs2", "rs3"],
        )

        assert outline_colors(fig) == [LEAD_OUTLINE] * 5 + [SECONDARY_OUTLINE] * 10

    def test_raises_for_invalid_highlight_snps(self, small_ld_matrix):
        """Test that invalid highlight_snps raises ValueError."""
        plotter = LDHeatmapPlotter()

        with pytest.raises(ValueError, match="highlight_snp"):
            plotter.plot_ld_heatmap(
                small_ld_matrix,
                highlight_snps=["invalid_snp"],
            )

    def test_metric_r2_label(self, small_ld_matrix):
        """Assert metric='r2' labels the colorbar R²."""
        plotter = LDHeatmapPlotter()
        fig = plotter.plot_ld_heatmap(small_ld_matrix, metric="r2")

        assert fig.get_axes()[1].get_ylabel() == "R²"

    def test_metric_dprime_label(self, small_ld_matrix):
        """Assert metric='dprime' labels the colorbar D'."""
        plotter = LDHeatmapPlotter()
        fig = plotter.plot_ld_heatmap(small_ld_matrix, metric="dprime")

        assert fig.get_axes()[1].get_ylabel() == "D'"

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

    def test_plotly_show_colorbar_false_hides_scale(self, small_ld_matrix):
        """Plotly honours show_colorbar=False rather than always showing it."""
        plotter = LDHeatmapPlotter(backend="plotly")
        fig = plotter.plot_ld_heatmap(small_ld_matrix, show_colorbar=False)

        assert [trace.showscale for trace in fig.data] == [False]

    def test_plotly_show_colorbar_true_shows_scale(self, small_ld_matrix):
        """Plotly still renders the scale when the caller asks for it."""
        plotter = LDHeatmapPlotter(backend="plotly")
        fig = plotter.plot_ld_heatmap(small_ld_matrix, show_colorbar=True)

        assert [trace.showscale for trace in fig.data] == [True]

    def test_plotly_colorbar_title_follows_metric(self, small_ld_matrix):
        """The label passed to add_colorbar reaches the Plotly colorbar title."""
        plotter = LDHeatmapPlotter(backend="plotly")

        r2_fig = plotter.plot_ld_heatmap(small_ld_matrix, metric="r2")
        dprime_fig = plotter.plot_ld_heatmap(small_ld_matrix, metric="dprime")

        assert r2_fig.data[0].colorbar.title.text == "R²"
        assert dprime_fig.data[0].colorbar.title.text == "D'"


class TestLDHeatmapBackends:
    """Tests for LDHeatmapPlotter with different backends."""

    @pytest.mark.parametrize("backend", BUILTIN_BACKENDS)
    def test_returns_the_backends_figure_type(self, backend, small_ld_matrix):
        """Each backend returns its own figure type."""
        plotter = LDHeatmapPlotter(backend=backend)
        fig = plotter.plot_ld_heatmap(small_ld_matrix)

        assert isinstance(fig, FIGURE_TYPES[backend])

    @pytest.mark.parametrize("backend", BUILTIN_BACKENDS)
    def test_nan_values_do_not_raise(self, backend, matrix_with_nan):
        """Each backend renders a matrix containing NaN."""
        plotter = LDHeatmapPlotter(backend=backend)
        plotter.plot_ld_heatmap(matrix_with_nan)


class TestLDHeatmapEdgeCases:
    """Tests for edge cases in LDHeatmapPlotter."""

    def test_single_snp_matrix_draws_one_cell(self):
        """Assert a 1x1 matrix draws its single cell and label."""
        plotter = LDHeatmapPlotter()
        matrix = pd.DataFrame([[1.0]], index=["rs1"], columns=["rs1"])

        fig = plotter.plot_ld_heatmap(matrix)

        np.testing.assert_array_equal(drawn_cells(fig), np.array([[1.0]]))
        assert tick_labels(fig) == (["rs1"], ["rs1"])

    def test_all_nan_matrix_draws_no_values(self):
        """Assert an all-NaN matrix leaves every drawn cell missing."""
        plotter = LDHeatmapPlotter()
        values = np.array(
            [
                [np.nan, np.nan],
                [np.nan, np.nan],
            ]
        )
        matrix = pd.DataFrame(values, index=["rs1", "rs2"], columns=["rs1", "rs2"])

        fig = plotter.plot_ld_heatmap(matrix)

        assert np.isnan(drawn_cells(fig)).all()

    def test_diagonal_only_values(self):
        """Assert an identity matrix draws ones on the diagonal and zeros below."""
        plotter = LDHeatmapPlotter()
        matrix = pd.DataFrame(
            np.eye(3),
            index=["rs1", "rs2", "rs3"],
            columns=["rs1", "rs2", "rs3"],
        )

        fig = plotter.plot_ld_heatmap(matrix)

        np.testing.assert_array_equal(
            drawn_cells(fig),
            np.array(
                [
                    [1.0, np.nan, np.nan],
                    [0.0, 1.0, np.nan],
                    [0.0, 0.0, 1.0],
                ]
            ),
        )

    def test_empty_highlight_snps_list(self, small_ld_matrix):
        """Assert an empty highlight_snps list outlines nothing."""
        plotter = LDHeatmapPlotter()
        fig = plotter.plot_ld_heatmap(small_ld_matrix, highlight_snps=[])

        assert outline_colors(fig) == []

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
        """Assert bokeh outlines the lead SNP's row and column in the lead colour."""
        from bokeh.models import GlyphRenderer, Rect

        plotter = LDHeatmapPlotter(backend="bokeh")
        fig = plotter.plot_ld_heatmap(sample_ld_matrix, lead_snp="rs3")

        outlines = [
            renderer.data_source.data
            for child in fig.children
            for renderer in child.renderers
            if isinstance(renderer, GlyphRenderer)
            and isinstance(renderer.glyph, Rect)
            and renderer.glyph.line_color == LEAD_SNP_HIGHLIGHT_COLOR
        ]
        centres = sorted((data["x"][0], data["y"][0]) for data in outlines)

        assert centres == [(0.0, 2.0), (1.0, 2.0), (2.0, 2.0), (2.0, 3.0), (2.0, 4.0)]

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
