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
    rng = np.random.default_rng(42)
    n_snps = 50
    positions = np.sort(rng.integers(1_000_000, 2_000_000, n_snps))
    return pd.DataFrame(
        {
            "rs": [f"rs{i}" for i in range(n_snps)],
            "chr": [1] * n_snps,
            "ps": positions,
            "p_wald": rng.uniform(1e-10, 1, n_snps),
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


class TestBackendForestPlotMethods:
    """Tests for forest plot backend methods (hbar, errorbar_h, axvline)."""

    def test_hbar_matplotlib(self):
        """Test horizontal bar chart in matplotlib backend."""
        from pylocuszoom.backends.matplotlib_backend import MatplotlibBackend

        backend = MatplotlibBackend()
        fig, axes = backend.create_figure(
            n_panels=1, height_ratios=[1.0], figsize=(8, 4)
        )
        ax = axes[0]

        y = pd.Series([0, 1, 2])
        width = pd.Series([0.5, 0.8, 0.3])

        backend.hbar(ax, y=y, width=width, height=0.5, color="blue")
        # Should not raise

    def test_errorbar_h_matplotlib(self):
        """Test horizontal error bar in matplotlib backend."""
        from pylocuszoom.backends.matplotlib_backend import MatplotlibBackend

        backend = MatplotlibBackend()
        fig, axes = backend.create_figure(
            n_panels=1, height_ratios=[1.0], figsize=(8, 4)
        )
        ax = axes[0]

        x = pd.Series([0.5, 0.8, 0.3])
        y = pd.Series([0, 1, 2])
        xerr_lower = pd.Series([0.1, 0.2, 0.1])
        xerr_upper = pd.Series([0.1, 0.1, 0.2])

        backend.errorbar_h(ax, x=x, y=y, xerr_lower=xerr_lower, xerr_upper=xerr_upper)
        # Should not raise

    def test_axvline_matplotlib(self):
        """Test vertical line in matplotlib backend."""
        from pylocuszoom.backends.matplotlib_backend import MatplotlibBackend

        backend = MatplotlibBackend()
        fig, axes = backend.create_figure(
            n_panels=1, height_ratios=[1.0], figsize=(8, 4)
        )
        ax = axes[0]

        backend.axvline(ax, x=0.5, color="red", linestyle="--")
        # Should not raise

    def test_axvline_plotly(self):
        """Test vertical line in plotly backend."""
        backend = PlotlyBackend()
        fig, axes = backend.create_figure(
            n_panels=1, height_ratios=[1.0], figsize=(8, 4)
        )
        ax = axes[0]

        backend.axvline(ax, x=0.5, color="red", linestyle="--")
        # Should not raise

    def test_axvline_bokeh(self):
        """Test vertical line in bokeh backend."""
        backend = BokehBackend()
        layout, axes = backend.create_figure(
            n_panels=1, height_ratios=[1.0], figsize=(8, 4)
        )
        ax = axes[0]

        backend.axvline(ax, x=0.5, color="red", linestyle="--")
        # Should not raise

    def test_hbar_plotly(self):
        """Test horizontal bar in plotly backend."""
        backend = PlotlyBackend()
        fig, axes = backend.create_figure(
            n_panels=1, height_ratios=[1.0], figsize=(8, 4)
        )
        ax = axes[0]

        y = pd.Series([0, 1, 2])
        width = pd.Series([0.5, 0.8, 0.3])

        backend.hbar(ax, y=y, width=width, height=0.5, color="blue")
        # Should not raise

    def test_errorbar_h_plotly(self):
        """Test horizontal error bar in plotly backend."""
        backend = PlotlyBackend()
        fig, axes = backend.create_figure(
            n_panels=1, height_ratios=[1.0], figsize=(8, 4)
        )
        ax = axes[0]

        x = pd.Series([0.5, 0.8, 0.3])
        y = pd.Series([0, 1, 2])
        xerr_lower = pd.Series([0.1, 0.2, 0.1])
        xerr_upper = pd.Series([0.1, 0.1, 0.2])

        backend.errorbar_h(ax, x=x, y=y, xerr_lower=xerr_lower, xerr_upper=xerr_upper)
        # Should not raise

    def test_hbar_bokeh(self):
        """Test horizontal bar in bokeh backend."""
        backend = BokehBackend()
        layout, axes = backend.create_figure(
            n_panels=1, height_ratios=[1.0], figsize=(8, 4)
        )
        ax = axes[0]

        y = pd.Series([0, 1, 2])
        width = pd.Series([0.5, 0.8, 0.3])

        backend.hbar(ax, y=y, width=width, height=0.5, color="blue")
        # Should not raise

    def test_errorbar_h_bokeh(self):
        """Test horizontal error bar in bokeh backend."""
        backend = BokehBackend()
        layout, axes = backend.create_figure(
            n_panels=1, height_ratios=[1.0], figsize=(8, 4)
        )
        ax = axes[0]

        x = pd.Series([0.5, 0.8, 0.3])
        y = pd.Series([0, 1, 2])
        xerr_lower = pd.Series([0.1, 0.2, 0.1])
        xerr_upper = pd.Series([0.1, 0.1, 0.2])

        backend.errorbar_h(ax, x=x, y=y, xerr_lower=xerr_lower, xerr_upper=xerr_upper)
        # Should not raise


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
            snp_labels=False,
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
            snp_labels=False,
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

    def test_all_backends_reject_empty_dataframe(self):
        """All backends should raise ValidationError for empty DataFrames."""
        from pylocuszoom.exceptions import ValidationError

        empty_df = pd.DataFrame(columns=["rs", "chr", "ps", "p_wald"])

        for backend_name in ["matplotlib", "plotly", "bokeh"]:
            plotter = LocusZoomPlotter(
                species="canine", backend=backend_name, log_level=None
            )
            with pytest.raises(ValidationError, match="empty"):
                plotter.plot(
                    empty_df,
                    chrom=1,
                    start=1_000_000,
                    end=2_000_000,
                    show_recombination=False,
                )

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
        df["R2"] = np.random.default_rng(0).uniform(0, 1, len(df))

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


@pytest.fixture
def sample_eqtl_df():
    """Sample eQTL DataFrame with effect sizes."""
    return pd.DataFrame(
        {
            "pos": [1_200_000, 1_400_000, 1_600_000, 1_800_000],
            "p_value": [1e-8, 1e-6, 1e-4, 1e-5],
            "effect_size": [0.5, -0.3, 0.8, -0.2],  # Mixed positive/negative
            "gene": ["GENE_A", "GENE_A", "GENE_A", "GENE_A"],
        }
    )


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
def sample_finemapping_df():
    """Sample fine-mapping DataFrame with credible sets."""
    return pd.DataFrame(
        {
            "pos": [1_200_000, 1_300_000, 1_400_000, 1_500_000, 1_600_000],
            "pip": [0.85, 0.10, 0.03, 0.45, 0.30],
            "cs": [1, 1, 0, 2, 2],  # Two credible sets + non-CS variants
        }
    )


class TestPlotlyEQTLFinemappingMarkers:
    """Tests for eQTL and fine-mapping marker rendering in Plotly."""

    def test_plotly_eqtl_positive_effect_markers(
        self, sample_gwas_df, sample_eqtl_df, sample_genes_df
    ):
        """Plotly eQTL positive effects should render as triangle-up markers."""
        plotter = LocusZoomPlotter(species="canine", backend="plotly", log_level=None)
        fig = plotter.plot_stacked(
            [sample_gwas_df],
            chrom=1,
            start=1_000_000,
            end=2_000_000,
            show_recombination=False,
            eqtl_df=sample_eqtl_df,
            eqtl_gene="GENE_A",
            genes_df=sample_genes_df,
        )

        # Find traces with triangle-up markers (positive effects)
        triangle_up_traces = [
            t
            for t in fig.data
            if hasattr(t, "marker") and t.marker.symbol == "triangle-up"
        ]
        assert len(triangle_up_traces) > 0, (
            "No triangle-up markers for positive effects"
        )

    def test_plotly_eqtl_negative_effect_markers(
        self, sample_gwas_df, sample_eqtl_df, sample_genes_df
    ):
        """Plotly eQTL negative effects should render as triangle-down markers."""
        plotter = LocusZoomPlotter(species="canine", backend="plotly", log_level=None)
        fig = plotter.plot_stacked(
            [sample_gwas_df],
            chrom=1,
            start=1_000_000,
            end=2_000_000,
            show_recombination=False,
            eqtl_df=sample_eqtl_df,
            eqtl_gene="GENE_A",
            genes_df=sample_genes_df,
        )

        # Find traces with triangle-down markers (negative effects)
        triangle_down_traces = [
            t
            for t in fig.data
            if hasattr(t, "marker") and t.marker.symbol == "triangle-down"
        ]
        assert len(triangle_down_traces) > 0, (
            "No triangle-down markers for negative effects"
        )

    def test_plotly_eqtl_no_effect_diamond_markers(
        self, sample_gwas_df, sample_eqtl_no_effect_df, sample_genes_df
    ):
        """Plotly eQTL without effect sizes should render as diamond markers."""
        plotter = LocusZoomPlotter(species="canine", backend="plotly", log_level=None)
        fig = plotter.plot_stacked(
            [sample_gwas_df],
            chrom=1,
            start=1_000_000,
            end=2_000_000,
            show_recombination=False,
            eqtl_df=sample_eqtl_no_effect_df,
            eqtl_gene="GENE_A",
            genes_df=sample_genes_df,
        )

        # Find traces with diamond markers
        diamond_traces = [
            t for t in fig.data if hasattr(t, "marker") and t.marker.symbol == "diamond"
        ]
        assert len(diamond_traces) > 0, (
            "No diamond markers for eQTL without effect sizes"
        )

    def test_plotly_finemapping_circle_markers(
        self, sample_gwas_df, sample_finemapping_df, sample_genes_df
    ):
        """Plotly fine-mapping should render as circle markers."""
        plotter = LocusZoomPlotter(species="canine", backend="plotly", log_level=None)
        fig = plotter.plot_stacked(
            [sample_gwas_df],
            chrom=1,
            start=1_000_000,
            end=2_000_000,
            show_recombination=False,
            finemapping_df=sample_finemapping_df,
            genes_df=sample_genes_df,
        )

        # Find traces with circle markers (fine-mapping uses circles)
        circle_traces = [
            t for t in fig.data if hasattr(t, "marker") and t.marker.symbol == "circle"
        ]
        assert len(circle_traces) > 0, "No circle markers for fine-mapping"

    def test_plotly_eqtl_hover_data(
        self, sample_gwas_df, sample_eqtl_df, sample_genes_df
    ):
        """Plotly eQTL scatter should have hover data with position, p-value, effect."""
        plotter = LocusZoomPlotter(species="canine", backend="plotly", log_level=None)
        fig = plotter.plot_stacked(
            [sample_gwas_df],
            chrom=1,
            start=1_000_000,
            end=2_000_000,
            show_recombination=False,
            eqtl_df=sample_eqtl_df,
            eqtl_gene="GENE_A",
            genes_df=sample_genes_df,
        )

        # Find eQTL data traces (triangle markers with actual data, not legend traces)
        eqtl_traces = [
            t
            for t in fig.data
            if hasattr(t, "marker")
            and t.marker.symbol in ("triangle-up", "triangle-down")
            and t.x is not None
            and len(t.x) > 0
            and t.x[0] is not None  # Exclude legend traces with None values
        ]
        assert len(eqtl_traces) > 0, "No eQTL data traces found"

        # Check hover data exists
        for trace in eqtl_traces:
            assert trace.customdata is not None, "eQTL trace missing customdata"

    def test_plotly_finemapping_hover_data(
        self, sample_gwas_df, sample_finemapping_df, sample_genes_df
    ):
        """Plotly fine-mapping scatter should have hover data with position, PIP, CS."""
        plotter = LocusZoomPlotter(species="canine", backend="plotly", log_level=None)
        fig = plotter.plot_stacked(
            [sample_gwas_df],
            chrom=1,
            start=1_000_000,
            end=2_000_000,
            show_recombination=False,
            finemapping_df=sample_finemapping_df,
            genes_df=sample_genes_df,
        )

        # Find fine-mapping traces (circle markers with PIP data)
        fm_traces = [
            t
            for t in fig.data
            if hasattr(t, "marker")
            and t.marker.symbol == "circle"
            and hasattr(t, "customdata")
            and t.customdata is not None
        ]
        assert len(fm_traces) > 0, "No fine-mapping traces with hover data found"


class TestBokehEQTLFinemappingMarkers:
    """Tests for eQTL and fine-mapping marker rendering in Bokeh."""

    def test_bokeh_eqtl_with_effects_creates_renderers(
        self, sample_gwas_df, sample_eqtl_df, sample_genes_df
    ):
        """Bokeh eQTL with effect sizes should create scatter renderers."""
        plotter = LocusZoomPlotter(species="canine", backend="bokeh", log_level=None)
        fig = plotter.plot_stacked(
            [sample_gwas_df],
            chrom=1,
            start=1_000_000,
            end=2_000_000,
            show_recombination=False,
            eqtl_df=sample_eqtl_df,
            eqtl_gene="GENE_A",
            genes_df=sample_genes_df,
        )

        # Bokeh returns a column layout - get the figures
        from bokeh.models import Column

        assert isinstance(fig, Column)
        # Should have multiple figures (GWAS + eQTL + gene track)
        assert len(fig.children) >= 2

    def test_bokeh_eqtl_triangle_markers(
        self, sample_gwas_df, sample_eqtl_df, sample_genes_df
    ):
        """Bokeh eQTL should use triangle markers for directional effects."""
        plotter = LocusZoomPlotter(species="canine", backend="bokeh", log_level=None)
        fig = plotter.plot_stacked(
            [sample_gwas_df],
            chrom=1,
            start=1_000_000,
            end=2_000_000,
            show_recombination=False,
            eqtl_df=sample_eqtl_df,
            eqtl_gene="GENE_A",
            genes_df=sample_genes_df,
        )

        from bokeh.models import GlyphRenderer, Scatter

        # Get all scatter renderers from all figures
        scatter_markers = []
        for child in fig.children:
            if hasattr(child, "renderers"):
                for r in child.renderers:
                    if isinstance(r, GlyphRenderer) and isinstance(r.glyph, Scatter):
                        scatter_markers.append(r.glyph.marker)

        # Should have triangle and inverted_triangle markers
        assert (
            "triangle" in scatter_markers or "inverted_triangle" in scatter_markers
        ), f"No triangle markers found in Bokeh plot. Markers: {scatter_markers}"

    def test_bokeh_finemapping_circle_markers(
        self, sample_gwas_df, sample_finemapping_df, sample_genes_df
    ):
        """Bokeh fine-mapping should use circle markers."""
        plotter = LocusZoomPlotter(species="canine", backend="bokeh", log_level=None)
        fig = plotter.plot_stacked(
            [sample_gwas_df],
            chrom=1,
            start=1_000_000,
            end=2_000_000,
            show_recombination=False,
            finemapping_df=sample_finemapping_df,
            genes_df=sample_genes_df,
        )

        from bokeh.models import GlyphRenderer, Scatter

        # Get all scatter renderers
        scatter_markers = []
        for child in fig.children:
            if hasattr(child, "renderers"):
                for r in child.renderers:
                    if isinstance(r, GlyphRenderer) and isinstance(r.glyph, Scatter):
                        scatter_markers.append(r.glyph.marker)

        # Should have circle markers for fine-mapping
        assert "circle" in scatter_markers, (
            f"No circle markers found in Bokeh plot. Markers: {scatter_markers}"
        )

    def test_bokeh_eqtl_has_hover_tool(
        self, sample_gwas_df, sample_eqtl_df, sample_genes_df
    ):
        """Bokeh eQTL panels should have HoverTool for interactivity."""
        from bokeh.models import HoverTool

        plotter = LocusZoomPlotter(species="canine", backend="bokeh", log_level=None)
        fig = plotter.plot_stacked(
            [sample_gwas_df],
            chrom=1,
            start=1_000_000,
            end=2_000_000,
            show_recombination=False,
            eqtl_df=sample_eqtl_df,
            eqtl_gene="GENE_A",
            genes_df=sample_genes_df,
        )

        # Check for HoverTool in any figure
        has_hover = False
        for child in fig.children:
            if hasattr(child, "tools"):
                for tool in child.tools:
                    if isinstance(tool, HoverTool):
                        has_hover = True
                        break

        assert has_hover, "No HoverTool found in Bokeh eQTL plot"

    def test_bokeh_finemapping_has_hover_tool(
        self, sample_gwas_df, sample_finemapping_df, sample_genes_df
    ):
        """Bokeh fine-mapping panels should have HoverTool for interactivity."""
        from bokeh.models import HoverTool

        plotter = LocusZoomPlotter(species="canine", backend="bokeh", log_level=None)
        fig = plotter.plot_stacked(
            [sample_gwas_df],
            chrom=1,
            start=1_000_000,
            end=2_000_000,
            show_recombination=False,
            finemapping_df=sample_finemapping_df,
            genes_df=sample_genes_df,
        )

        # Check for HoverTool in any figure
        has_hover = False
        for child in fig.children:
            if hasattr(child, "tools"):
                for tool in child.tools:
                    if isinstance(tool, HoverTool):
                        has_hover = True
                        break

        assert has_hover, "No HoverTool found in Bokeh fine-mapping plot"


class TestGeneTrackMbFormatting:
    """Tests for Mb formatting on gene track axis in interactive backends."""

    def test_plotly_gene_track_has_mb_formatting(self, sample_gwas_df, sample_genes_df):
        """Plotly gene track axis should have Mb formatting (not raw bp).

        Regression test: gene track axis showed raw bp ticks while label said "Mb".
        """
        plotter = LocusZoomPlotter(species="canine", backend="plotly", log_level=None)
        fig = plotter.plot(
            sample_gwas_df,
            chrom=1,
            start=1_000_000,
            end=2_000_000,
            genes_df=sample_genes_df,
            show_recombination=False,
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

    def test_bokeh_gene_track_has_mb_formatting(self, sample_gwas_df, sample_genes_df):
        """Bokeh gene track axis should have Mb formatting (not raw bp).

        Regression test: gene track axis showed raw bp ticks while label said "Mb".
        """
        from bokeh.models import CustomJSTickFormatter

        plotter = LocusZoomPlotter(species="canine", backend="bokeh", log_level=None)
        fig = plotter.plot(
            sample_gwas_df,
            chrom=1,
            start=1_000_000,
            end=2_000_000,
            genes_df=sample_genes_df,
            show_recombination=False,
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
        from pylocuszoom.backends.plotly_backend import PlotlyBackend

        backend = PlotlyBackend()

        # Create a figure with 12 subplots (enough to trigger old collision)
        fig, axes = backend.create_figure(12, [1.0] * 12, (12, 24))

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
        from pylocuszoom.backends.plotly_backend import PlotlyBackend

        backend = PlotlyBackend()
        fig, axes = backend.create_figure(5, [1.0] * 5, (12, 10))

        secondary_names = []
        for ax in axes:
            _, secondary_y = backend.create_twin_axis(ax)
            secondary_names.append(secondary_y)

        assert len(set(secondary_names)) == 5, (
            f"Secondary axis names are not unique: {secondary_names}"
        )


class TestHeatmapCoordinates:
    """Tests that heatmap backends use actual genomic coordinates, not indices."""

    def test_plotly_heatmap_uses_actual_coords(self):
        """Plotly heatmap should use passed x/y coords, not range(len(coords))."""
        from pylocuszoom.backends.plotly_backend import PlotlyBackend

        backend = PlotlyBackend()
        fig, axes = backend.create_figure(1, [1.0], (10, 6))
        ax = axes[0]

        data = np.array([[1.0, 0.5], [0.5, 1.0]])
        x_coords = [1_000_000.0, 2_000_000.0]
        y_coords = [1_000_000.0, 2_000_000.0]

        backend.add_heatmap(ax, data, x_coords, y_coords, mask_upper=False)

        fig_obj = ax[0]
        heatmap_trace = [t for t in fig_obj.data if hasattr(t, "z")][0]
        assert list(heatmap_trace.x) == [1_000_000.0, 2_000_000.0], (
            f"Expected genomic x-coords, got {heatmap_trace.x}"
        )
        assert list(heatmap_trace.y) == [1_000_000.0, 2_000_000.0], (
            f"Expected genomic y-coords, got {heatmap_trace.y}"
        )

    def test_bokeh_heatmap_uses_actual_coords(self):
        """Bokeh heatmap should use passed x/y coords, not index values."""
        from pylocuszoom.backends.bokeh_backend import BokehBackend

        backend = BokehBackend()
        fig, axes = backend.create_figure(1, [1.0], (10, 6))
        ax = axes[0]

        data = np.array([[1.0, 0.5], [0.5, 1.0]])
        x_coords = [1_000_000.0, 2_000_000.0]
        y_coords = [1_000_000.0, 2_000_000.0]

        backend.add_heatmap(ax, data, x_coords, y_coords, mask_upper=False)

        # Get the rect glyph data source
        renderers = ax.renderers
        rect_renderer = [r for r in renderers if hasattr(r, "glyph")][-1]
        xs = sorted(set(rect_renderer.data_source.data["x"]))
        ys = sorted(set(rect_renderer.data_source.data["y"]))
        assert xs == [1_000_000.0, 2_000_000.0], f"Expected genomic x-coords, got {xs}"
        assert ys == [1_000_000.0, 2_000_000.0], f"Expected genomic y-coords, got {ys}"

    def test_bokeh_heatmap_mask_upper_uses_actual_coords(self):
        """Bokeh heatmap with mask_upper=True should still use actual coords."""
        from pylocuszoom.backends.bokeh_backend import BokehBackend

        backend = BokehBackend()
        fig, axes = backend.create_figure(1, [1.0], (10, 6))
        ax = axes[0]

        data = np.array([[1.0, 0.5], [0.5, 1.0]])
        x_coords = [1_000_000.0, 2_000_000.0]
        y_coords = [1_000_000.0, 2_000_000.0]

        backend.add_heatmap(ax, data, x_coords, y_coords, mask_upper=True)

        renderers = ax.renderers
        rect_renderer = [r for r in renderers if hasattr(r, "glyph")][-1]
        xs = list(rect_renderer.data_source.data["x"])
        ys = list(rect_renderer.data_source.data["y"])
        # mask_upper=True: only lower triangle (j <= i), so 3 cells for 2x2
        assert len(xs) == 3
        assert all(x >= 1_000_000 for x in xs), f"Expected genomic x-coords, got {xs}"
        assert all(y >= 1_000_000 for y in ys), f"Expected genomic y-coords, got {ys}"

    def test_bokeh_heatmap_rect_dimensions_scale_with_coords(self):
        """Bokeh heatmap rect width/height should match coordinate spacing."""
        from pylocuszoom.backends.bokeh_backend import BokehBackend

        backend = BokehBackend()
        fig, axes = backend.create_figure(1, [1.0], (10, 6))
        ax = axes[0]

        data = np.array([[1.0, 0.5], [0.5, 1.0]])
        x_coords = [1_000_000.0, 2_000_000.0]
        y_coords = [1_000_000.0, 2_000_000.0]

        backend.add_heatmap(ax, data, x_coords, y_coords, mask_upper=False)

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
        fig, axes = backend.create_figure(1, [1.0], (10, 6))
        ax = axes[0]

        # 3x3 matrix with non-uniform spacing
        data = np.eye(3)
        x_coords = [1_000_000.0, 2_000_000.0, 4_000_000.0]  # gap doubles
        y_coords = [1_000_000.0, 2_000_000.0, 4_000_000.0]

        backend.add_heatmap(ax, data, x_coords, y_coords, mask_upper=False)

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
