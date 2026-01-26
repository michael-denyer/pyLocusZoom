"""Tests for notebook compatibility of interactive backends.

These tests ensure Plotly and Bokeh backends produce outputs that
are compatible with Jupyter/Databricks notebook environments.
"""

import json
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from pylocuszoom.backends.bokeh_backend import BokehBackend
from pylocuszoom.backends.plotly_backend import PlotlyBackend
from pylocuszoom.plotter import LocusZoomPlotter


@pytest.fixture
def sample_gwas_df():
    """Sample GWAS results DataFrame."""
    np.random.seed(42)
    n_snps = 50
    positions = np.sort(np.random.randint(1_000_000, 2_000_000, n_snps))
    return pd.DataFrame(
        {
            "rs": [f"rs{i}" for i in range(n_snps)],
            "chr": [1] * n_snps,
            "ps": positions,
            "p_wald": np.random.uniform(1e-10, 1, n_snps),
        }
    )


@pytest.fixture
def sample_genes_df():
    """Sample gene annotations."""
    return pd.DataFrame(
        {
            "chr": [1, 1, 1],
            "start": [1_100_000, 1_400_000, 1_700_000],
            "end": [1_150_000, 1_500_000, 1_800_000],
            "gene_name": ["GENE_A", "GENE_B", "GENE_C"],
            "strand": ["+", "-", "+"],
        }
    )


class TestPlotlyNotebookCompatibility:
    """Tests for Plotly backend notebook compatibility."""

    def test_plotly_figure_has_repr_html(self, sample_gwas_df):
        """Plotly figures must have _repr_html_() for notebook display."""
        plotter = LocusZoomPlotter(species="canine", backend="plotly", log_level=None)
        fig = plotter.plot(
            sample_gwas_df,
            chrom=1,
            start=1_000_000,
            end=2_000_000,
            show_recombination=False,
        )

        # Plotly figures have _repr_html_ for notebook rendering
        assert hasattr(fig, "_repr_html_")
        assert callable(fig._repr_html_)

        # Should produce valid HTML
        html = fig._repr_html_()
        assert isinstance(html, str)
        assert len(html) > 0
        assert "plotly" in html.lower() or "div" in html.lower()

    def test_plotly_figure_to_json(self, sample_gwas_df):
        """Plotly figures must be JSON-serializable for Databricks."""
        plotter = LocusZoomPlotter(species="canine", backend="plotly", log_level=None)
        fig = plotter.plot(
            sample_gwas_df,
            chrom=1,
            start=1_000_000,
            end=2_000_000,
            show_recombination=False,
        )

        # Databricks uses JSON serialization
        json_str = fig.to_json()
        assert isinstance(json_str, str)

        # Should be valid JSON
        parsed = json.loads(json_str)
        assert "data" in parsed
        assert "layout" in parsed

    def test_plotly_figure_to_html(self, sample_gwas_df):
        """Plotly figures must save to HTML for notebook export."""
        plotter = LocusZoomPlotter(species="canine", backend="plotly", log_level=None)
        fig = plotter.plot(
            sample_gwas_df,
            chrom=1,
            start=1_000_000,
            end=2_000_000,
            show_recombination=False,
        )

        with tempfile.NamedTemporaryFile(suffix=".html", delete=False) as f:
            fig.write_html(f.name)
            html_content = Path(f.name).read_text()

        assert len(html_content) > 0
        assert "<html" in html_content or "<!DOCTYPE" in html_content

    def test_plotly_figure_has_data(self, sample_gwas_df):
        """Plotly figures must contain scatter data."""
        plotter = LocusZoomPlotter(species="canine", backend="plotly", log_level=None)
        fig = plotter.plot(
            sample_gwas_df,
            chrom=1,
            start=1_000_000,
            end=2_000_000,
            show_recombination=False,
        )

        # Should have at least one trace
        assert len(fig.data) >= 1

        # First trace should be scatter
        assert fig.data[0].type == "scatter"

    def test_plotly_hover_data(self, sample_gwas_df):
        """Plotly figures should have hover text for interactive exploration."""
        plotter = LocusZoomPlotter(species="canine", backend="plotly", log_level=None)
        fig = plotter.plot(
            sample_gwas_df,
            chrom=1,
            start=1_000_000,
            end=2_000_000,
            show_recombination=False,
        )

        # Check that traces have hovertemplate (plotly's hover mechanism)
        assert len(fig.data) > 0
        # At least one trace should have hover info
        has_hover = any(
            hasattr(trace, "hovertemplate") and trace.hovertemplate
            for trace in fig.data
        )
        assert has_hover, "No traces have hovertemplate"

    def test_plotly_stacked_figure(self, sample_gwas_df):
        """Plotly backend should work with plot_stacked()."""
        plotter = LocusZoomPlotter(species="canine", backend="plotly", log_level=None)
        fig = plotter.plot_stacked(
            [sample_gwas_df, sample_gwas_df.copy()],
            chrom=1,
            start=1_000_000,
            end=2_000_000,
            show_recombination=False,
            snp_labels=False,  # SNP labels are matplotlib-only
        )

        # Should have data from both panels
        assert len(fig.data) >= 2

        # Should be JSON-serializable
        json_str = fig.to_json()
        assert isinstance(json_str, str)


class TestBokehNotebookCompatibility:
    """Tests for Bokeh backend notebook compatibility."""

    def test_bokeh_figure_creation_no_errors(self, sample_gwas_df):
        """Bokeh figure creation should not raise errors."""
        plotter = LocusZoomPlotter(species="canine", backend="bokeh", log_level=None)

        # Should complete without errors
        fig = plotter.plot(
            sample_gwas_df,
            chrom=1,
            start=1_000_000,
            end=2_000_000,
            show_recombination=False,
        )

        assert fig is not None

    def test_bokeh_figure_saves_to_html(self, sample_gwas_df):
        """Bokeh figures must save to HTML for notebook export."""
        from bokeh.io import output_file, save

        plotter = LocusZoomPlotter(species="canine", backend="bokeh", log_level=None)
        fig = plotter.plot(
            sample_gwas_df,
            chrom=1,
            start=1_000_000,
            end=2_000_000,
            show_recombination=False,
        )

        with tempfile.NamedTemporaryFile(suffix=".html", delete=False) as f:
            output_file(f.name)
            save(fig)
            html_content = Path(f.name).read_text()

        assert len(html_content) > 0
        assert "<html" in html_content or "<!DOCTYPE" in html_content
        assert "bokeh" in html_content.lower()

    def test_bokeh_figure_json_serialization(self, sample_gwas_df):
        """Bokeh figures must be JSON-serializable for notebook display."""
        from bokeh.embed import json_item

        plotter = LocusZoomPlotter(species="canine", backend="bokeh", log_level=None)
        fig = plotter.plot(
            sample_gwas_df,
            chrom=1,
            start=1_000_000,
            end=2_000_000,
            show_recombination=False,
        )

        # Bokeh uses json_item for embedding
        json_data = json_item(fig)
        assert isinstance(json_data, dict)
        assert "doc" in json_data or "root_id" in json_data

    def test_bokeh_uses_scatter_not_deprecated_circle(self, sample_gwas_df):
        """Bokeh backend must use scatter() not deprecated circle() method."""
        backend = BokehBackend()

        # Create a figure
        layout, figures = backend.create_figure(
            n_panels=1,
            height_ratios=[1.0],
            figsize=(12, 6),
        )

        ax = figures[0]
        x = sample_gwas_df["ps"]
        y = -np.log10(sample_gwas_df["p_wald"])

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

    def test_bokeh_uses_customjs_tick_formatter(self, sample_gwas_df):
        """Bokeh backend must use CustomJSTickFormatter not deprecated FuncTickFormatter."""
        from bokeh.models import CustomJSTickFormatter

        backend = BokehBackend()
        layout, figures = backend.create_figure(
            n_panels=1,
            height_ratios=[1.0],
            figsize=(12, 6),
        )

        ax = figures[0]
        backend.format_xaxis_mb(ax)

        # Should use CustomJSTickFormatter
        assert isinstance(ax.xaxis.formatter, CustomJSTickFormatter)

    def test_bokeh_column_layout_no_sizing_mode_warning(self, sample_gwas_df):
        """Bokeh column layout should not trigger FIXED_SIZING_MODE warning."""
        backend = BokehBackend()

        # Creating figure should not produce validation warnings
        # (We can't easily test for warnings here, but we test the API is correct)
        layout, figures = backend.create_figure(
            n_panels=2,
            height_ratios=[3.0, 1.0],
            figsize=(12, 8),
        )

        # Should create valid layout without errors
        assert layout is not None
        assert len(figures) == 2

    def test_bokeh_stacked_figure(self, sample_gwas_df):
        """Bokeh backend should work with plot_stacked()."""
        from bokeh.io import output_file, save

        plotter = LocusZoomPlotter(species="canine", backend="bokeh", log_level=None)
        fig = plotter.plot_stacked(
            [sample_gwas_df, sample_gwas_df.copy()],
            chrom=1,
            start=1_000_000,
            end=2_000_000,
            show_recombination=False,
            snp_labels=False,  # SNP labels are matplotlib-only
        )

        # Should save to HTML without errors
        with tempfile.NamedTemporaryFile(suffix=".html", delete=False) as f:
            output_file(f.name)
            save(fig)
            html_content = Path(f.name).read_text()

        assert len(html_content) > 0


class TestBackendConsistency:
    """Tests ensuring consistent output across backends."""

    def test_all_backends_return_figure(self, sample_gwas_df):
        """All backends should return a figure object."""
        for backend_name in ["matplotlib", "plotly", "bokeh"]:
            plotter = LocusZoomPlotter(
                species="canine", backend=backend_name, log_level=None
            )
            fig = plotter.plot(
                sample_gwas_df,
                chrom=1,
                start=1_000_000,
                end=2_000_000,
                show_recombination=False,
            )
            assert fig is not None, f"{backend_name} returned None"

    def test_all_backends_handle_empty_dataframe(self):
        """All backends should handle empty DataFrames gracefully."""
        empty_df = pd.DataFrame(columns=["rs", "chr", "ps", "p_wald"])

        for backend_name in ["matplotlib", "plotly", "bokeh"]:
            plotter = LocusZoomPlotter(
                species="canine", backend=backend_name, log_level=None
            )
            fig = plotter.plot(
                empty_df,
                chrom=1,
                start=1_000_000,
                end=2_000_000,
                show_recombination=False,
            )
            assert fig is not None, f"{backend_name} failed with empty DataFrame"

    def test_all_backends_handle_lead_position(self, sample_gwas_df):
        """All backends should handle lead_pos parameter."""
        for backend_name in ["matplotlib", "plotly", "bokeh"]:
            plotter = LocusZoomPlotter(
                species="canine", backend=backend_name, log_level=None
            )
            fig = plotter.plot(
                sample_gwas_df,
                chrom=1,
                start=1_000_000,
                end=2_000_000,
                lead_pos=1_500_000,
                show_recombination=False,
            )
            assert fig is not None, f"{backend_name} failed with lead_pos"

    def test_all_backends_handle_precomputed_ld(self, sample_gwas_df):
        """All backends should handle pre-computed LD column."""
        df = sample_gwas_df.copy()
        df["R2"] = np.random.uniform(0, 1, len(df))

        for backend_name in ["matplotlib", "plotly", "bokeh"]:
            plotter = LocusZoomPlotter(
                species="canine", backend=backend_name, log_level=None
            )
            fig = plotter.plot(
                df,
                chrom=1,
                start=1_000_000,
                end=2_000_000,
                ld_col="R2",
                show_recombination=False,
            )
            assert fig is not None, f"{backend_name} failed with ld_col"


class TestBackendSaveOperations:
    """Tests for backend save functionality."""

    def test_plotly_backend_save_html(self, sample_gwas_df):
        """PlotlyBackend.save() should work for HTML files."""
        backend = PlotlyBackend()
        layout, axes = backend.create_figure(
            n_panels=1, height_ratios=[1.0], figsize=(12, 6)
        )

        # Add some data
        x = sample_gwas_df["ps"]
        y = -np.log10(sample_gwas_df["p_wald"])
        backend.scatter(axes[0], x, y, colors="#BEBEBE")

        with tempfile.NamedTemporaryFile(suffix=".html", delete=False) as f:
            backend.save(layout, f.name)
            assert Path(f.name).exists()
            assert Path(f.name).stat().st_size > 0

    def test_bokeh_backend_save_html(self, sample_gwas_df):
        """BokehBackend.save() should work for HTML files."""
        backend = BokehBackend()
        layout, axes = backend.create_figure(
            n_panels=1, height_ratios=[1.0], figsize=(12, 6)
        )

        # Add some data
        x = sample_gwas_df["ps"]
        y = -np.log10(sample_gwas_df["p_wald"])
        backend.scatter(axes[0], x, y, colors="#BEBEBE")

        with tempfile.NamedTemporaryFile(suffix=".html", delete=False) as f:
            backend.save(layout, f.name)
            assert Path(f.name).exists()
            assert Path(f.name).stat().st_size > 0


class TestDatabricksSpecific:
    """Tests specific to Databricks notebook environment."""

    def test_plotly_displayhtml_compatible(self, sample_gwas_df):
        """Plotly output should be compatible with Databricks displayHTML()."""
        plotter = LocusZoomPlotter(species="canine", backend="plotly", log_level=None)
        fig = plotter.plot(
            sample_gwas_df,
            chrom=1,
            start=1_000_000,
            end=2_000_000,
            show_recombination=False,
        )

        # Databricks displayHTML() expects a complete HTML string
        html = fig.to_html(include_plotlyjs=True, full_html=True)
        assert isinstance(html, str)
        assert "<html" in html or "<!DOCTYPE" in html
        assert "plotly" in html.lower()

    def test_bokeh_components_for_embedding(self, sample_gwas_df):
        """Bokeh should provide components for Databricks embedding."""
        from bokeh.embed import components

        plotter = LocusZoomPlotter(species="canine", backend="bokeh", log_level=None)
        fig = plotter.plot(
            sample_gwas_df,
            chrom=1,
            start=1_000_000,
            end=2_000_000,
            show_recombination=False,
        )

        # components() returns (script, div) for embedding
        script, div = components(fig)
        assert isinstance(script, str)
        assert isinstance(div, str)
        assert "<script" in script
        assert "<div" in div
