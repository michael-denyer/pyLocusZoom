"""Tests for regional plots rendered through each backend."""

import pandas as pd
import pytest
from hypothesis import given
from hypothesis import settings as hyp_settings

from pylocuszoom import DisplayConfig, PanelInputs
from pylocuszoom.backends import BUILTIN_BACKENDS
from pylocuszoom.plotter import LocusZoomPlotter
from tests.conftest import FIGURE_TYPES
from tests.strategies import gwas_dataframes

PANEL_COUNTS = {
    "matplotlib": lambda fig: len(fig.get_axes()),
    "plotly": lambda fig: sum(
        key.startswith("yaxis") for key in fig.layout.to_plotly_json()
    ),
    "bokeh": lambda fig: len(fig.children),
}
"""How many stacked panels a figure carries, in each backend's own terms."""


class TestBackendIntegration:
    """Tests for backend protocol integration."""

    def test_default_backend_draws_a_matplotlib_figure(self, tiny_regional_gwas_df):
        """An unconfigured plotter renders through matplotlib."""
        plotter = LocusZoomPlotter(species="canine")

        fig = plotter.plot(
            tiny_regional_gwas_df,
            chrom=1,
            start=1000000,
            end=2000000,
            display=DisplayConfig(show_recombination=False),
        )

        assert isinstance(fig, FIGURE_TYPES["matplotlib"])

    def test_plotly_backend_creates_figure(self, tiny_regional_gwas_df):
        """plot() with backend='plotly' produces a plotly Figure.

        This also implicitly confirms plot() routes through the backend
        protocol: if plot() bypassed the backend and called matplotlib
        directly, the returned object would be a matplotlib Figure and
        this isinstance check would fail.
        """
        import plotly.graph_objects as go

        plotter = LocusZoomPlotter(species="canine", backend="plotly")

        fig = plotter.plot(
            tiny_regional_gwas_df,
            chrom=1,
            start=1000000,
            end=2000000,
            display=DisplayConfig(show_recombination=False),
        )

        assert isinstance(fig, go.Figure)

    def test_matplotlib_plot_renders_expected_artists(self, tiny_regional_gwas_df):
        """plot() renders scatter points, a significance line, and axis labels.

        Assertions query the rendered axes directly (scatter collections,
        line objects, y-label text) rather than mocking backend methods
        and counting calls — per CLAUDE.md's observable-outputs rule.
        """
        plotter = LocusZoomPlotter(species="canine")

        fig = plotter.plot(
            tiny_regional_gwas_df,
            chrom=1,
            start=1000000,
            end=2000000,
            display=DisplayConfig(show_recombination=False),
        )
        ax = fig.axes[0]

        # At least one scatter collection (association points)
        assert len(ax.collections) >= 1, (
            "expected at least one scatter collection on association axes"
        )

        # At least one line (the significance threshold)
        lines = ax.get_lines()
        assert len(lines) >= 1, "expected a significance threshold line"

        # Y-axis labelled with -log10(p)
        ylabel = ax.get_ylabel()
        assert "log" in ylabel.lower() or "-log" in ylabel.lower(), (
            f"expected -log10(p) y-label, got {ylabel!r}"
        )

        # X-limits cover the requested region
        xlim = ax.get_xlim()
        assert xlim[0] <= 1_000_000 and xlim[1] >= 2_000_000

    def test_plot_stacked_renders_two_panels(self, tiny_regional_gwas_df):
        """plot_stacked() with two GWAS inputs produces two association panels."""
        plotter = LocusZoomPlotter(species="canine")

        fig = plotter.plot_stacked(
            [tiny_regional_gwas_df, tiny_regional_gwas_df.copy()],
            chrom=1,
            start=1000000,
            end=2000000,
            display=DisplayConfig(show_recombination=False),
        )
        # Each GWAS gets its own association axes; gene track axes
        # may also be present but at minimum we need two scatter axes.
        scatter_axes = [ax for ax in fig.axes if ax.collections]
        assert len(scatter_axes) >= 2, (
            f"expected >=2 axes with scatter data, got {len(scatter_axes)}"
        )


class TestBackendEQTLFinemapping:
    """Tests for eQTL and fine-mapping support across all backends."""

    @pytest.fixture
    def sample_eqtl_df_no_effect(self):
        """Sample eQTL DataFrame without effect sizes."""
        return pd.DataFrame(
            {
                "pos": [1200000, 1400000, 1600000],
                "p_value": [1e-6, 1e-4, 0.01],
                "gene": ["GENE1", "GENE1", "GENE1"],
            }
        )

    @pytest.mark.parametrize("backend", BUILTIN_BACKENDS)
    @pytest.mark.parametrize(
        "eqtl_fixture,with_finemapping,expected_panels",
        [
            pytest.param("sample_eqtl_df", False, 2, id="eqtl-with-effects"),
            pytest.param("sample_eqtl_df_no_effect", False, 2, id="eqtl-no-effects"),
            pytest.param(None, True, 2, id="finemapping"),
            pytest.param("sample_eqtl_df", True, 3, id="eqtl-and-finemapping"),
        ],
    )
    def test_optional_panels_render_on_every_backend(
        self,
        request,
        backend,
        eqtl_fixture,
        with_finemapping,
        expected_panels,
        small_regional_gwas_df,
        sample_finemapping_df,
    ):
        """Every backend draws the association panel plus each optional panel asked for."""
        plotter = LocusZoomPlotter(species=None, backend=backend, log_level=None)

        panels = {}
        if eqtl_fixture is not None:
            panels["eqtl_df"] = request.getfixturevalue(eqtl_fixture)
            panels["eqtl_gene"] = "GENE1"
        if with_finemapping:
            panels["finemapping_df"] = sample_finemapping_df

        fig = plotter.plot_stacked(
            [small_regional_gwas_df],
            chrom=1,
            start=1000000,
            end=2000000,
            panels=PanelInputs(**panels),
            display=DisplayConfig(show_recombination=False),
        )

        assert isinstance(fig, FIGURE_TYPES[backend])
        assert PANEL_COUNTS[backend](fig) == expected_panels

    def test_plot_accepts_eqtl_and_finemapping_panels(
        self, small_regional_gwas_df, sample_eqtl_df, sample_finemapping_df
    ):
        """plot() carries the same optional panels as plot_stacked()."""
        plotter = LocusZoomPlotter(species=None, backend="matplotlib", log_level=None)

        fig = plotter.plot(
            small_regional_gwas_df,
            chrom=1,
            start=1000000,
            end=2000000,
            display=DisplayConfig(show_recombination=False),
            panels=PanelInputs(
                eqtl_df=sample_eqtl_df,
                eqtl_gene="GENE1",
                finemapping_df=sample_finemapping_df,
            ),
        )

        axes = fig.get_axes()
        assert len(axes) == 3
        assert axes[1].get_ylabel() == "PIP"
        assert "eQTL" in axes[2].get_ylabel()

    def test_eqtl_chr_filtering(self, small_regional_gwas_df):
        """Drop eQTLs on another chromosome even when their position is in range."""
        plotter = LocusZoomPlotter(species=None, backend="matplotlib", log_level=None)

        eqtl_df = pd.DataFrame(
            {
                "pos": [1200000, 1400000, 1600000],  # All in region 1e6-2e6
                "p_value": [1e-6, 1e-4, 0.01],
                "gene": ["GENE1", "GENE1", "GENE1"],
                "effect_size": [0.5, -0.3, 0.8],
                "chr": ["1", "2", "1"],  # Middle one is on chr2
            }
        )

        fig = plotter.plot_stacked(
            [small_regional_gwas_df],
            chrom=1,
            start=1000000,
            end=2000000,
            display=DisplayConfig(show_recombination=False),
            panels=PanelInputs(eqtl_df=eqtl_df, eqtl_gene="GENE1"),
        )

        eqtl_ax = fig.get_axes()[1]
        assert "eQTL" in eqtl_ax.get_ylabel()
        drawn = {
            float(x)
            for collection in eqtl_ax.collections
            for x, _ in collection.get_offsets()
        }
        assert drawn == {1200000.0, 1600000.0}

    def test_eqtl_gene_without_gene_column_raises(self, small_regional_gwas_df):
        """eqtl_gene on a frame with no gene column is an error, not a warning.

        Drawing every eQTL unfiltered after the caller asked for one gene is
        the silent-failure shape the library rejects.
        """
        from pylocuszoom.eqtl import EQTLValidationError

        plotter = LocusZoomPlotter(species=None, backend="matplotlib", log_level=None)
        eqtl_df_no_gene_col = pd.DataFrame(
            {
                "pos": [1200000, 1400000, 1600000],
                "p_value": [1e-6, 1e-4, 0.01],
            }
        )

        with pytest.raises(EQTLValidationError, match="gene"):
            plotter.plot_stacked(
                [small_regional_gwas_df],
                chrom=1,
                start=1000000,
                end=2000000,
                display=DisplayConfig(show_recombination=False),
                panels=PanelInputs(eqtl_df=eqtl_df_no_gene_col, eqtl_gene="GENE1"),
            )

    def test_eqtl_zero_pvalue_is_dropped(self, small_regional_gwas_df):
        """A zero eQTL p-value is outside the strict (0, 1] domain and not drawn."""
        plotter = LocusZoomPlotter(species=None, backend="matplotlib", log_level=None)
        eqtl_df = pd.DataFrame(
            {
                "pos": [1200000, 1400000, 1600000],
                "p_value": [1e-6, 0.0, 0.01],
            }
        )

        fig = plotter.plot_stacked(
            [small_regional_gwas_df],
            chrom=1,
            start=1000000,
            end=2000000,
            display=DisplayConfig(show_recombination=False),
            panels=PanelInputs(eqtl_df=eqtl_df),
        )

        eqtl_ax = fig.get_axes()[1]
        drawn = sum(len(c.get_offsets()) for c in eqtl_ax.collections)
        assert drawn == 2


class TestPlotterProperties:
    """Property-based tests for plotter crash-resistance."""

    @pytest.mark.parametrize("backend", BUILTIN_BACKENDS)
    @hyp_settings(max_examples=10, deadline=None)
    @given(gwas_dataframes(min_snps=3, max_snps=30))
    def test_plot_renders_valid_data(self, backend, df):
        """plot() renders any valid GWAS data on any backend without crashing."""
        plotter = LocusZoomPlotter(species="canine", backend=backend)

        fig = plotter.plot(
            df,
            chrom=df["chr"].iloc[0],
            start=int(df["pos"].min()),
            end=int(df["pos"].max()),
            display=DisplayConfig(show_recombination=False),
        )

        assert isinstance(fig, FIGURE_TYPES[backend])
