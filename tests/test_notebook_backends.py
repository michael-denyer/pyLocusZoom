"""Tests for notebook compatibility of interactive backends.

These tests ensure Plotly and Bokeh backends produce outputs that
are compatible with Jupyter/Databricks notebook environments.
"""

import json

import numpy as np
import pandas as pd
import pytest

from pylocuszoom import DisplayConfig, LDConfig, PanelInputs
from pylocuszoom.backends import BUILTIN_BACKENDS, get_backend
from pylocuszoom.backends.bokeh_backend import BokehBackend
from pylocuszoom.backends.composition import lower_triangle
from pylocuszoom.backends.plotly_backend import PlotlyBackend
from pylocuszoom.colors import LD_HEATMAP_COLORS
from pylocuszoom.plotter import LocusZoomPlotter
from tests.conftest import FIGURE_TYPES
from tests.figure_probes import INTERACTIVE_BACKENDS, PROBES


class TestBackendForestPlotMethods:
    """Tests for forest plot backend methods (errorbar_h, axvline)."""

    @staticmethod
    def _draw_errorbar_h(backend, ax):
        backend.errorbar_h(
            ax,
            x=pd.Series([0.5, 0.8, 0.3]),
            y=pd.Series([0, 1, 2]),
            xerr_lower=pd.Series([0.1, 0.2, 0.1]),
            xerr_upper=pd.Series([0.1, 0.1, 0.2]),
        )

    @staticmethod
    def _draw_axvline(backend, ax):
        backend.axvline(ax, x=0.5, color="red", linestyle="--")

    @pytest.mark.parametrize("backend_name", BUILTIN_BACKENDS)
    @pytest.mark.parametrize("primitive", ["_draw_errorbar_h", "_draw_axvline"])
    def test_forest_primitive_does_not_raise(self, backend_name, primitive):
        """Every backend accepts the two primitives a forest plot draws."""
        backend = get_backend(backend_name)
        _, axes = backend.create_figure(height_ratios=[1.0], figsize=(8, 4))

        getattr(self, primitive)(backend, axes[0])


@pytest.fixture
def regional_figure(backend_name, regional_gwas_df):
    """The single-panel regional figure the export tests inspect."""
    plotter = LocusZoomPlotter(species="canine", backend=backend_name, log_level=None)
    return plotter.plot(
        regional_gwas_df,
        chrom=1,
        start=1_000_000,
        end=2_000_000,
        display=DisplayConfig(show_recombination=False),
    )


@pytest.mark.parametrize("backend_name", INTERACTIVE_BACKENDS)
class TestNotebookExport:
    """Every interactive backend produces what a notebook needs to display it."""

    def test_figure_exports_a_standalone_html_document(
        self, backend_name, regional_figure
    ):
        """Databricks displayHTML() is handed one string carrying its own library."""
        html = PROBES[backend_name].standalone_html(regional_figure)

        assert "<html" in html or "<!DOCTYPE" in html
        assert backend_name in html.lower()

    def test_figure_serialises_to_json(self, backend_name, regional_figure):
        """The front end is handed a payload it can post as JSON."""
        payload = PROBES[backend_name].json_payload(regional_figure)

        assert isinstance(payload, dict)
        assert payload
        json.dumps(payload)

    def test_stacked_figure_exports_the_same_way(self, backend_name, regional_gwas_df):
        """Two stacked panels export as one document, not two."""
        plotter = LocusZoomPlotter(
            species="canine", backend=backend_name, log_level=None
        )
        fig = plotter.plot_stacked(
            [regional_gwas_df, regional_gwas_df.copy()],
            chrom=1,
            start=1_000_000,
            end=2_000_000,
            display=DisplayConfig(show_recombination=False, snp_labels=False),
        )

        html = PROBES[backend_name].standalone_html(fig)

        assert "<html" in html or "<!DOCTYPE" in html

    def test_figure_carries_hover(self, backend_name, regional_figure):
        """A regional plot is explorable, not a static image."""
        assert PROBES[backend_name].has_hover(regional_figure)


class TestPlotlyNotebookCompatibility:
    """Plotly-only notebook affordances."""

    def test_plotly_figure_has_repr_html(self, regional_gwas_df):
        """Jupyter renders a plotly figure through _repr_html_()."""
        plotter = LocusZoomPlotter(species="canine", backend="plotly", log_level=None)
        fig = plotter.plot(
            regional_gwas_df,
            chrom=1,
            start=1_000_000,
            end=2_000_000,
            display=DisplayConfig(show_recombination=False),
        )

        html = fig._repr_html_()

        assert isinstance(html, str)
        assert "plotly" in html.lower() or "div" in html.lower()

    def test_plotly_figure_has_data(self, regional_gwas_df):
        """Plotly figures must contain scatter data."""
        plotter = LocusZoomPlotter(species="canine", backend="plotly", log_level=None)
        fig = plotter.plot(
            regional_gwas_df,
            chrom=1,
            start=1_000_000,
            end=2_000_000,
            display=DisplayConfig(show_recombination=False),
        )

        # Should have at least one trace
        assert len(fig.data) >= 1

        # First trace should be scatter
        assert fig.data[0].type == "scatter"


class TestBokehNotebookCompatibility:
    """Bokeh-only notebook affordances."""

    def test_bokeh_components_for_embedding(self, regional_gwas_df):
        """Bokeh should provide components for Databricks embedding."""
        from bokeh.embed import components

        plotter = LocusZoomPlotter(species="canine", backend="bokeh", log_level=None)
        fig = plotter.plot(
            regional_gwas_df,
            chrom=1,
            start=1_000_000,
            end=2_000_000,
            display=DisplayConfig(show_recombination=False),
        )

        script, div = components(fig)

        assert "<script" in script
        assert "<div" in div

    def test_bokeh_uses_scatter_not_deprecated_circle(self, regional_gwas_df):
        """Bokeh backend must use scatter() not deprecated circle() method."""
        backend = BokehBackend()

        # Create a figure
        layout, figures = backend.create_figure(
            height_ratios=[1.0],
            figsize=(12, 6),
        )

        ax = figures[0]
        x = regional_gwas_df["pos"]
        y = -np.log10(regional_gwas_df["p_value"])

        # scatter() should work without deprecation warning
        backend.scatter(
            ax=ax,
            x=x,
            y=y,
            colors="#BEBEBE",
            sizes=60,
        )

        # Should have renderers
        assert len(ax.renderers) > 0

    def test_bokeh_uses_customjs_tick_formatter(self):
        """Bokeh backend must use CustomJSTickFormatter not deprecated FuncTickFormatter."""
        from bokeh.models import CustomJSTickFormatter

        backend = BokehBackend()
        layout, figures = backend.create_figure(
            height_ratios=[1.0],
            figsize=(12, 6),
        )

        ax = figures[0]
        backend.format_xaxis_mb(ax)

        # Should use CustomJSTickFormatter
        assert isinstance(ax.xaxis.formatter, CustomJSTickFormatter)

    def test_bokeh_column_layout_no_sizing_mode_warning(self):
        """A stacked layout passes bokeh's own integrity check.

        Bokeh's FIXED_SIZING_MODE warning fires when a layout asks for fixed
        sizing without a width and height, which is what a notebook user sees
        as a blank or clipped column. ``check_integrity`` is the same pass
        bokeh runs on save, so it fails here rather than in the browser.
        """
        from bokeh.core.validation import check_integrity

        backend = BokehBackend()
        layout, figures = backend.create_figure(
            height_ratios=[3.0, 1.0],
            figsize=(12, 8),
        )

        issues = check_integrity([layout])

        assert [w.name for w in issues.warning] == []
        assert [e.name for e in issues.error] == []
        assert len(figures) == 2


class TestBackendConsistency:
    """Every built-in backend renders the same regional plot request."""

    @pytest.mark.parametrize("backend_name", BUILTIN_BACKENDS)
    def test_backend_returns_its_own_figure_type(self, backend_name, regional_gwas_df):
        """Each backend returns the figure type its library defines."""
        plotter = LocusZoomPlotter(
            species="canine", backend=backend_name, log_level=None
        )
        fig = plotter.plot(
            regional_gwas_df,
            chrom=1,
            start=1_000_000,
            end=2_000_000,
            display=DisplayConfig(show_recombination=False),
        )
        assert isinstance(fig, FIGURE_TYPES[backend_name])

    @pytest.mark.parametrize("backend_name", BUILTIN_BACKENDS)
    def test_backend_rejects_empty_dataframe(self, backend_name):
        """Each backend raises ValidationError for an empty frame."""
        from pylocuszoom.exceptions import ValidationError

        empty_df = pd.DataFrame(columns=["rs", "chr", "pos", "p_value"])
        plotter = LocusZoomPlotter(
            species="canine", backend=backend_name, log_level=None
        )
        with pytest.raises(ValidationError, match="empty"):
            plotter.plot(
                empty_df,
                chrom=1,
                start=1_000_000,
                end=2_000_000,
                display=DisplayConfig(show_recombination=False),
            )

    @pytest.mark.parametrize("backend_name", BUILTIN_BACKENDS)
    def test_backend_handles_lead_position(self, backend_name, regional_gwas_df):
        """Each backend renders a figure when lead_pos is given."""
        plotter = LocusZoomPlotter(
            species="canine", backend=backend_name, log_level=None
        )
        fig = plotter.plot(
            regional_gwas_df,
            chrom=1,
            start=1_000_000,
            end=2_000_000,
            display=DisplayConfig(show_recombination=False),
            ld=LDConfig(lead_pos=1_500_000),
        )
        assert isinstance(fig, FIGURE_TYPES[backend_name])

    @pytest.mark.parametrize("backend_name", BUILTIN_BACKENDS)
    def test_backend_handles_precomputed_ld(self, backend_name, regional_gwas_df):
        """Each backend renders a figure from a precomputed LD column."""
        df = regional_gwas_df.assign(
            R2=np.random.default_rng(0).uniform(0, 1, len(regional_gwas_df))
        )
        plotter = LocusZoomPlotter(
            species="canine", backend=backend_name, log_level=None
        )
        fig = plotter.plot(
            df,
            chrom=1,
            start=1_000_000,
            end=2_000_000,
            display=DisplayConfig(show_recombination=False),
            ld=LDConfig(ld_col="R2"),
        )
        assert isinstance(fig, FIGURE_TYPES[backend_name])


@pytest.fixture
def sample_eqtl_no_effect_df():
    """Sample eQTL DataFrame without effect sizes."""
    return pd.DataFrame(
        {
            "pos": [1_200_000, 1_400_000, 1_600_000],
            "p_value": [1e-8, 1e-6, 1e-4],
            "gene": ["GENE_A", "GENE_A", "GENE_A"],
        }
    )


@pytest.fixture
def eqtl_figure(backend_name, regional_gwas_df, sample_eqtl_df, sample_genes_df):
    """GWAS, an eQTL panel with signed effects, and a gene track."""
    return _stacked(
        backend_name,
        regional_gwas_df,
        sample_genes_df,
        eqtl_df=sample_eqtl_df,
        eqtl_gene="GENE_A",
    )


@pytest.fixture
def finemapping_figure(
    backend_name, regional_gwas_df, sample_finemapping_df, sample_genes_df
):
    """GWAS, a fine-mapping panel and a gene track."""
    return _stacked(
        backend_name,
        regional_gwas_df,
        sample_genes_df,
        finemapping_df=sample_finemapping_df,
    )


def _stacked(backend_name, gwas_df, genes_df, **panels):
    plotter = LocusZoomPlotter(species="canine", backend=backend_name, log_level=None)
    return plotter.plot_stacked(
        [gwas_df],
        chrom=1,
        start=1_000_000,
        end=2_000_000,
        display=DisplayConfig(show_recombination=False),
        panels=PanelInputs(genes_df=genes_df, **panels),
    )


@pytest.mark.parametrize("backend_name", INTERACTIVE_BACKENDS)
class TestOptionalPanelMarkers:
    """Marker choice is a library-independent decision, so both backends make it."""

    def test_effect_direction_gets_its_own_triangle(self, backend_name, eqtl_figure):
        """An eQTL panel points a triangle up or down per the sign of the effect.

        Both directions are asserted. The bokeh half of this used to pass with
        the negative-effect glyph never drawn.
        """
        symbols = PROBES[backend_name].marker_symbols(eqtl_figure)

        assert "triangle-up" in symbols
        assert "triangle-down" in symbols

    def test_eqtl_without_effects_draws_no_triangles(
        self, backend_name, regional_gwas_df, sample_eqtl_no_effect_df, sample_genes_df
    ):
        """With no effect column there is no direction, so the marker is neutral."""
        fig = _stacked(
            backend_name,
            regional_gwas_df,
            sample_genes_df,
            eqtl_df=sample_eqtl_no_effect_df,
            eqtl_gene="GENE_A",
        )

        symbols = PROBES[backend_name].marker_symbols(fig)

        assert "diamond" in symbols
        assert "triangle-up" not in symbols
        assert "triangle-down" not in symbols

    def test_finemapping_uses_circles(self, backend_name, finemapping_figure):
        """A PIP has no direction, so fine-mapping draws circles."""
        assert "circle" in PROBES[backend_name].marker_symbols(finemapping_figure)

    def test_eqtl_panel_carries_hover(self, backend_name, eqtl_figure):
        """An eQTL point is worth hovering: it carries a gene, a p-value, an effect."""
        assert PROBES[backend_name].has_hover(eqtl_figure)

    def test_finemapping_panel_carries_hover(self, backend_name, finemapping_figure):
        """A fine-mapping point carries its PIP and credible set."""
        assert PROBES[backend_name].has_hover(finemapping_figure)


class TestGeneTrackMbFormatting:
    """Tests for Mb formatting on gene track axis in interactive backends."""

    def test_plotly_gene_track_has_mb_formatting(
        self, regional_gwas_df, sample_genes_df
    ):
        """Plotly gene track axis should have Mb formatting (not raw bp).

        Regression test: gene track axis showed raw bp ticks while label said "Mb".
        """
        plotter = LocusZoomPlotter(species="canine", backend="plotly", log_level=None)
        fig = plotter.plot(
            regional_gwas_df,
            chrom=1,
            start=1_000_000,
            end=2_000_000,
            display=DisplayConfig(show_recombination=False),
            panels=PanelInputs(genes_df=sample_genes_df),
        )

        # With gene track, row 2 is the gene track axis, so xaxis2 is its x-axis.
        gene_track_xaxis = fig.layout.xaxis2
        assert gene_track_xaxis.ticksuffix == " Mb", (
            f"Gene track axis ticksuffix is {gene_track_xaxis.ticksuffix!r}, "
            "expected ' Mb'"
        )
        assert gene_track_xaxis.ticktext, "Gene track axis has no Mb tick labels"
        assert gene_track_xaxis.ticktext[0] == "1.00", (
            f"Gene track ticks start at {gene_track_xaxis.ticktext[0]!r}, "
            "expected '1.00' for a region starting at 1 Mb"
        )

    def test_bokeh_gene_track_has_mb_formatting(
        self, regional_gwas_df, sample_genes_df
    ):
        """Bokeh gene track axis should have Mb formatting (not raw bp).

        Regression test: gene track axis showed raw bp ticks while label said "Mb".
        """
        from bokeh.models import CustomJSTickFormatter

        plotter = LocusZoomPlotter(species="canine", backend="bokeh", log_level=None)
        fig = plotter.plot(
            regional_gwas_df,
            chrom=1,
            start=1_000_000,
            end=2_000_000,
            display=DisplayConfig(show_recombination=False),
            panels=PanelInputs(genes_df=sample_genes_df),
        )

        # Bokeh layout contains multiple figures - find the gene track figure
        # The gene track is typically the second figure (index 1)
        gene_track_fig = fig.children[1] if len(fig.children) > 1 else fig.children[0]

        # Check that the x-axis has CustomJSTickFormatter for Mb formatting
        assert isinstance(gene_track_fig.xaxis.formatter, CustomJSTickFormatter), (
            f"Gene track x-axis formatter is {type(gene_track_fig.xaxis.formatter)}, "
            "expected CustomJSTickFormatter for Mb formatting"
        )


class TestPlotlySecondaryAxisNaming:
    """Tests for Plotly secondary axis naming to avoid collisions."""

    def test_plotly_secondary_axis_offset_avoids_primary_collision(self):
        """Secondary axis suffix should not collide with primary subplot axes.

        Regression test: create_twin_axis used offset of 10, which collides
        with primary yaxis10+ when figures have 10+ subplots.
        """

        backend = PlotlyBackend()

        # Create a figure with 12 subplots (enough to trigger old collision)
        fig, axes = backend.create_figure([1.0] * 12, (12, 24))

        # Create twin axis for the first subplot
        _, secondary_y = backend.create_twin_axis(axes[0])

        # The secondary axis name should use high offset (100+), not 10+
        # y100 should not collide with any primary axis (y1 through y12)
        suffix = int(secondary_y.replace("y", ""))
        assert suffix >= 100, (
            f"Secondary axis suffix {suffix} is too low, "
            "will collide with primary axes in large subplot figures"
        )

    def test_plotly_secondary_axes_unique_across_subplots(self):
        """Each subplot's secondary axis should have a unique name."""

        backend = PlotlyBackend()
        fig, axes = backend.create_figure([1.0] * 5, (12, 10))

        secondary_names = []
        for ax in axes:
            _, secondary_y = backend.create_twin_axis(ax)
            secondary_names.append(secondary_y)

        assert len(set(secondary_names)) == 5, (
            f"Secondary axis names are not unique: {secondary_names}"
        )


class TestHeatmapCoordinates:
    """Tests that heatmap backends use actual genomic coordinates, not indices."""

    @pytest.mark.parametrize("backend_name", INTERACTIVE_BACKENDS)
    def test_heatmap_uses_the_coordinates_it_was_given(self, backend_name):
        """A heatmap is drawn at genomic positions, not at matrix indices."""
        backend = get_backend(backend_name)
        _, axes = backend.create_figure([1.0], (10, 6))
        ax = axes[0]
        coords = [1_000_000.0, 2_000_000.0]

        backend.add_heatmap(
            ax,
            np.array([[1.0, 0.5], [0.5, 1.0]]),
            coords,
            coords,
            cmap_colors=LD_HEATMAP_COLORS,
        )

        assert PROBES[backend_name].heatmap_coords(ax) == (coords, coords)

    def test_bokeh_heatmap_lower_triangle_uses_actual_coords(self):
        """Bokeh heatmap of a lower-triangle matrix should still use actual coords."""
        from pylocuszoom.backends.bokeh_backend import BokehBackend

        backend = BokehBackend()
        fig, axes = backend.create_figure([1.0], (10, 6))
        ax = axes[0]

        data = np.array([[1.0, 0.5], [0.5, 1.0]])
        x_coords = [1_000_000.0, 2_000_000.0]
        y_coords = [1_000_000.0, 2_000_000.0]

        backend.add_heatmap(
            ax, lower_triangle(data), x_coords, y_coords, cmap_colors=LD_HEATMAP_COLORS
        )

        renderers = ax.renderers
        rect_renderer = [r for r in renderers if hasattr(r, "glyph")][-1]
        xs = list(rect_renderer.data_source.data["x"])
        ys = list(rect_renderer.data_source.data["y"])
        # Lower triangle (j <= i) is 3 cells for a 2x2 matrix.
        assert len(xs) == 3
        assert all(x >= 1_000_000 for x in xs), f"Expected genomic x-coords, got {xs}"
        assert all(y >= 1_000_000 for y in ys), f"Expected genomic y-coords, got {ys}"

    def test_bokeh_heatmap_rect_dimensions_scale_with_coords(self):
        """Bokeh heatmap rect width/height should match coordinate spacing."""
        from pylocuszoom.backends.bokeh_backend import BokehBackend

        backend = BokehBackend()
        fig, axes = backend.create_figure([1.0], (10, 6))
        ax = axes[0]

        data = np.array([[1.0, 0.5], [0.5, 1.0]])
        x_coords = [1_000_000.0, 2_000_000.0]
        y_coords = [1_000_000.0, 2_000_000.0]

        backend.add_heatmap(ax, data, x_coords, y_coords, cmap_colors=LD_HEATMAP_COLORS)

        renderers = ax.renderers
        rect_renderer = [r for r in renderers if hasattr(r, "glyph")][-1]
        glyph = rect_renderer.glyph
        # Width/height are now per-cell column references ("w", "h")
        assert glyph.width == "w", f"Expected column ref 'w', got {glyph.width}"
        assert glyph.height == "h", f"Expected column ref 'h', got {glyph.height}"
        # Verify the data source contains correct cell sizes
        source = rect_renderer.data_source
        widths = source.data["w"]
        heights = source.data["h"]
        assert all(w == 1_000_000.0 for w in widths), (
            f"Expected widths 1M, got {widths}"
        )
        assert all(h == 1_000_000.0 for h in heights), (
            f"Expected heights 1M, got {heights}"
        )

    def test_bokeh_heatmap_nonuniform_spacing(self):
        """Bokeh heatmap per-cell sizing handles non-uniform coordinate spacing."""
        from pylocuszoom.backends.bokeh_backend import BokehBackend

        backend = BokehBackend()
        fig, axes = backend.create_figure([1.0], (10, 6))
        ax = axes[0]

        # 3x3 matrix with non-uniform spacing
        data = np.eye(3)
        x_coords = [1_000_000.0, 2_000_000.0, 4_000_000.0]  # gap doubles
        y_coords = [1_000_000.0, 2_000_000.0, 4_000_000.0]

        backend.add_heatmap(ax, data, x_coords, y_coords, cmap_colors=LD_HEATMAP_COLORS)

        renderers = ax.renderers
        rect_renderer = [r for r in renderers if hasattr(r, "glyph")][-1]
        source = rect_renderer.data_source
        widths = list(source.data["w"])

        # First cell: midpoint boundary at 1.5M, extrapolated left to 0.5M → width 1M
        # Second cell: boundaries at 1.5M and 3M → width 1.5M
        # Third cell: boundary at 3M, extrapolated right to 5M → width 2M
        # Pattern repeats for each row (9 cells total)
        assert len(widths) == 9
        # Widths should NOT all be equal (non-uniform spacing)
        assert len(set(round(w, 1) for w in widths)) > 1, (
            f"Expected non-uniform widths, got {widths}"
        )
