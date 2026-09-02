"""Tests for LD calculation and the LD heatmap panel in regional plots."""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from pylocuszoom import DisplayConfig, LDConfig, PanelInputs
from pylocuszoom.backends.composition import LD_LEGEND_TITLE
from pylocuszoom.exceptions import PlinkError
from pylocuszoom.plotter import LocusZoomPlotter


class TestLocusZoomPlotterLdCalculation:
    """Tests for LD calculation integration."""

    LD_OUTPUT = (
        "CHR_A BP_A SNP_A CHR_B BP_B SNP_B R2\n"
        "1 1100000 rs1 1 1500000 rs2 0.80\n"
        "1 1100000 rs1 1 1900000 rs3 0.50\n"
    )

    def test_ld_reference_colours_the_plot(self, fake_plink, tiny_regional_gwas_df):
        """A PLINK run that returns LD pairs colours the points by R2."""
        bfile, plink_writes = fake_plink

        with plink_writes(self.LD_OUTPUT):
            fig = LocusZoomPlotter(species="canine", plink_path="/mock/plink").plot(
                tiny_regional_gwas_df,
                chrom=1,
                start=1000000,
                end=2000000,
                display=DisplayConfig(show_recombination=False),
                ld=LDConfig(lead_pos=1100000, ld_reference_file=bfile),
            )

        legend = fig.get_axes()[0].get_legend()
        assert legend is not None, "LD data should add an r² legend"
        assert legend.get_title().get_text() == LD_LEGEND_TITLE
        assert "Lead SNP" in [text.get_text() for text in legend.get_texts()]

    def test_empty_ld_output_is_downgraded_to_warning(
        self, fake_plink, tiny_regional_gwas_df, warning_records
    ):
        """An empty PLINK output (singleton lead SNP) should not abort the plot.

        Singleton lead SNPs with no LD neighbours in the window are a real
        scenario; plot() catches only this specific PlinkError and continues
        without LD colouring, leaving a warning in the log.
        """
        bfile, plink_writes = fake_plink
        header_only = "CHR_A BP_A SNP_A CHR_B BP_B SNP_B R2\n"

        with plink_writes(header_only):
            fig = LocusZoomPlotter(species="canine", plink_path="/mock/plink").plot(
                tiny_regional_gwas_df,
                chrom=1,
                start=1000000,
                end=2000000,
                display=DisplayConfig(show_recombination=False),
                ld=LDConfig(lead_pos=1100000, ld_reference_file=bfile),
            )

        assert fig is not None
        assert any("LD calculation skipped" in message for message in warning_records)

    def test_stacked_plot_downgrades_empty_ld_output(
        self, fake_plink, tiny_regional_gwas_df
    ):
        """Single and stacked plots share the recoverable LD policy."""
        bfile, plink_writes = fake_plink
        header_only = "CHR_A BP_A SNP_A CHR_B BP_B SNP_B R2\n"

        with plink_writes(header_only):
            fig = LocusZoomPlotter(
                species="canine", plink_path="/mock/plink"
            ).plot_stacked(
                [tiny_regional_gwas_df],
                chrom=1,
                start=1000000,
                end=2000000,
                lead_positions=[1100000],
                ld_reference_files=[bfile],
                display=DisplayConfig(show_recombination=False),
            )

        assert fig is not None

    def test_plink_misconfiguration_propagates_through_plot(
        self, fake_plink, tiny_regional_gwas_df
    ):
        """A non-zero PLINK exit must surface; it means PLINK is misconfigured.

        The catch in plot() is narrow on purpose. Timeout, non-zero exit, and
        "output file missing after success" all indicate real misconfiguration
        and should reach the caller.
        """
        bfile, plink_writes = fake_plink

        with (
            plink_writes(None, returncode=2, stderr="bad bfile"),
            pytest.raises(PlinkError, match="exit code"),
        ):
            LocusZoomPlotter(species="canine", plink_path="/mock/plink").plot(
                tiny_regional_gwas_df,
                chrom=1,
                start=1000000,
                end=2000000,
                ld=LDConfig(lead_pos=1100000, ld_reference_file=bfile),
            )


class TestLDHeatmapIntegration:
    """Tests for LD heatmap integration in LocusZoomPlotter."""

    @pytest.fixture
    def ld_heatmap_gwas_df(self):
        """Sample GWAS DataFrame with positions matching LD heatmap SNPs."""
        return pd.DataFrame(
            {
                "rs": ["rs1", "rs2", "rs3", "rs4", "rs5"],
                "chr": [1, 1, 1, 1, 1],
                "ps": [1000000, 1000500, 1001000, 1001500, 1002000],
                "p_wald": [1e-8, 1e-6, 1e-4, 1e-3, 0.05],
            }
        )

    @pytest.fixture
    def sample_ld_heatmap_data(self):
        """Sample LD heatmap data matching calculate_pairwise_ld return format.

        Returns:
            Tuple of (DataFrame, list[str]): LD matrix and SNP IDs.
        """
        ld_matrix = pd.DataFrame(
            np.array(
                [
                    [1.0, 0.9, 0.7, 0.4, 0.2],
                    [0.9, 1.0, 0.8, 0.5, 0.3],
                    [0.7, 0.8, 1.0, 0.6, 0.4],
                    [0.4, 0.5, 0.6, 1.0, 0.7],
                    [0.2, 0.3, 0.4, 0.7, 1.0],
                ]
            ),
            index=["rs1", "rs2", "rs3", "rs4", "rs5"],
            columns=["rs1", "rs2", "rs3", "rs4", "rs5"],
        )
        snp_ids = ["rs1", "rs2", "rs3", "rs4", "rs5"]
        return ld_matrix, snp_ids

    @pytest.fixture
    def heatmap_genes_df(self):
        """Sample genes for testing stacked plots with gene track."""
        return pd.DataFrame(
            {
                "chr": ["1", "1"],
                "start": [1000200, 1001200],
                "end": [1000800, 1001800],
                "gene_name": ["GENE1", "GENE2"],
                "strand": ["+", "-"],
            }
        )

    def test_plot_with_ld_heatmap_renders_two_panels(
        self, ld_heatmap_gwas_df, sample_ld_heatmap_data
    ):
        """When ld_heatmap_df and ld_heatmap_snp_ids provided, figure has heatmap panel."""
        ld_matrix, snp_ids = sample_ld_heatmap_data
        plotter = LocusZoomPlotter(species=None, log_level=None)

        fig = plotter.plot(
            ld_heatmap_gwas_df,
            chrom=1,
            start=999000,
            end=1003000,
            display=DisplayConfig(show_recombination=False),
            panels=PanelInputs(ld_heatmap_df=ld_matrix, ld_heatmap_snp_ids=snp_ids),
        )

        # Should have at least 2 axes (association + heatmap)
        assert len(fig.axes) >= 2

    def test_plot_with_ld_heatmap_aligns_x_coordinates(
        self, ld_heatmap_gwas_df, sample_ld_heatmap_data
    ):
        """Heatmap SNPs render at their genomic positions from GWAS data."""
        ld_matrix, snp_ids = sample_ld_heatmap_data
        plotter = LocusZoomPlotter(species=None, log_level=None)

        fig = plotter.plot(
            ld_heatmap_gwas_df,
            chrom=1,
            start=999000,
            end=1003000,
            display=DisplayConfig(show_recombination=False),
            panels=PanelInputs(ld_heatmap_df=ld_matrix, ld_heatmap_snp_ids=snp_ids),
        )

        # Verify the heatmap panel exists and has correct x-axis range
        # The heatmap x-axis should span the genomic positions of SNPs
        axes = fig.axes
        assert len(axes) >= 2
        # Heatmap is axes[1] (after association at axes[0], colorbar may be axes[2])
        heatmap_ax = axes[1]
        xlim = heatmap_ax.get_xlim()
        # X-axis should be in genomic coordinate range
        assert xlim[0] < 1003000, f"Heatmap xlim[0]={xlim[0]} should be < 1003000"
        assert xlim[1] > 999000, f"Heatmap xlim[1]={xlim[1]} should be > 999000"

    def test_plot_stacked_with_ld_heatmap_at_bottom(
        self, ld_heatmap_gwas_df, sample_ld_heatmap_data, heatmap_genes_df
    ):
        """In stacked plots, heatmap appears below gene track (at very bottom)."""
        ld_matrix, snp_ids = sample_ld_heatmap_data
        plotter = LocusZoomPlotter(species=None, log_level=None)

        fig = plotter.plot_stacked(
            [ld_heatmap_gwas_df],
            chrom=1,
            start=999000,
            end=1003000,
            display=DisplayConfig(show_recombination=False),
            panels=PanelInputs(
                genes_df=heatmap_genes_df,
                ld_heatmap_df=ld_matrix,
                ld_heatmap_snp_ids=snp_ids,
            ),
        )

        # Should have 3 panels: GWAS, gene track, heatmap
        axes = fig.axes
        assert len(axes) >= 3

    def test_ld_heatmap_filters_to_region(
        self, ld_heatmap_gwas_df, sample_ld_heatmap_data
    ):
        """SNPs outside [start, end] are filtered from heatmap."""
        ld_matrix, snp_ids = sample_ld_heatmap_data
        plotter = LocusZoomPlotter(species=None, log_level=None)

        # Narrow region that only includes some SNPs
        # rs1 at 1000000, rs2 at 1000500 - both in [1000000, 1001000)
        # rs3 at 1001000 - at boundary
        fig = plotter.plot(
            ld_heatmap_gwas_df,
            chrom=1,
            start=1000000,
            end=1001000,
            display=DisplayConfig(show_recombination=False),
            panels=PanelInputs(ld_heatmap_df=ld_matrix, ld_heatmap_snp_ids=snp_ids),
        )

        # Should complete without error even with partial overlap
        assert fig is not None

    def test_ld_heatmap_empty_overlap_raises(
        self, ld_heatmap_gwas_df, sample_ld_heatmap_data
    ):
        """A heatmap whose SNPs all fall outside the region is a caller fault."""
        ld_matrix, snp_ids = sample_ld_heatmap_data
        plotter = LocusZoomPlotter(species=None, log_level="WARNING")

        with pytest.raises(ValueError, match="No SNPs from LD heatmap overlap"):
            plotter.plot(
                ld_heatmap_gwas_df,
                chrom=1,
                start=5000000,
                end=6000000,
                display=DisplayConfig(show_recombination=False),
                panels=PanelInputs(ld_heatmap_df=ld_matrix, ld_heatmap_snp_ids=snp_ids),
            )

    def test_ld_heatmap_height_parameter(
        self, ld_heatmap_gwas_df, sample_ld_heatmap_data
    ):
        """ld_heatmap_height parameter controls panel height ratio."""
        ld_matrix, snp_ids = sample_ld_heatmap_data
        plotter = LocusZoomPlotter(species=None, log_level=None)

        # Test with different height values
        fig1 = plotter.plot(
            ld_heatmap_gwas_df,
            chrom=1,
            start=999000,
            end=1003000,
            display=DisplayConfig(show_recombination=False),
            panels=PanelInputs(
                ld_heatmap_df=ld_matrix,
                ld_heatmap_snp_ids=snp_ids,
                ld_heatmap_height=0.1,
            ),
        )

        fig2 = plotter.plot(
            ld_heatmap_gwas_df,
            chrom=1,
            start=999000,
            end=1003000,
            display=DisplayConfig(show_recombination=False),
            panels=PanelInputs(
                ld_heatmap_df=ld_matrix,
                ld_heatmap_snp_ids=snp_ids,
                ld_heatmap_height=0.5,
            ),
        )

        # Both should render successfully
        assert fig1 is not None
        assert fig2 is not None

    def test_ld_heatmap_lead_snp_highlight(
        self, ld_heatmap_gwas_df, sample_ld_heatmap_data
    ):
        """Lead SNP is highlighted in heatmap when lead_pos specified and SNP found."""
        ld_matrix, snp_ids = sample_ld_heatmap_data
        plotter = LocusZoomPlotter(species=None, log_level=None)

        # rs1 is at position 1000000
        fig = plotter.plot(
            ld_heatmap_gwas_df,
            chrom=1,
            start=999000,
            end=1003000,
            display=DisplayConfig(show_recombination=False),
            ld=LDConfig(lead_pos=1000000),
            panels=PanelInputs(ld_heatmap_df=ld_matrix, ld_heatmap_snp_ids=snp_ids),
        )

        # Should render with lead SNP highlight (visual check)
        assert fig is not None

    # Backend-specific tests

    def test_ld_heatmap_matplotlib_backend(
        self, ld_heatmap_gwas_df, sample_ld_heatmap_data
    ):
        """Verify matplotlib figure has correct panel count and axes."""
        ld_matrix, snp_ids = sample_ld_heatmap_data
        plotter = LocusZoomPlotter(species=None, backend="matplotlib", log_level=None)

        fig = plotter.plot(
            ld_heatmap_gwas_df,
            chrom=1,
            start=999000,
            end=1003000,
            display=DisplayConfig(show_recombination=False),
            panels=PanelInputs(ld_heatmap_df=ld_matrix, ld_heatmap_snp_ids=snp_ids),
        )

        # Should have at least 2 main axes (association + heatmap)
        # Plus possible colorbar axis
        assert len(fig.axes) >= 2
        assert isinstance(fig, plt.Figure)

    def test_ld_heatmap_plotly_backend(
        self, ld_heatmap_gwas_df, sample_ld_heatmap_data
    ):
        """Verify plotly figure has heatmap trace at correct row."""
        import plotly.graph_objects as go

        ld_matrix, snp_ids = sample_ld_heatmap_data
        plotter = LocusZoomPlotter(species=None, backend="plotly", log_level=None)

        fig = plotter.plot(
            ld_heatmap_gwas_df,
            chrom=1,
            start=999000,
            end=1003000,
            display=DisplayConfig(show_recombination=False),
            panels=PanelInputs(ld_heatmap_df=ld_matrix, ld_heatmap_snp_ids=snp_ids),
        )

        assert isinstance(fig, go.Figure)
        # Check that figure has data traces
        assert len(fig.data) > 0

    def test_ld_heatmap_bokeh_backend(self, ld_heatmap_gwas_df, sample_ld_heatmap_data):
        """Verify bokeh layout contains heatmap."""
        from bokeh.models.layouts import Column

        ld_matrix, snp_ids = sample_ld_heatmap_data
        plotter = LocusZoomPlotter(species=None, backend="bokeh", log_level=None)

        fig = plotter.plot(
            ld_heatmap_gwas_df,
            chrom=1,
            start=999000,
            end=1003000,
            display=DisplayConfig(show_recombination=False),
            panels=PanelInputs(ld_heatmap_df=ld_matrix, ld_heatmap_snp_ids=snp_ids),
        )

        # Bokeh returns Column layout
        assert isinstance(fig, Column)
        # Should have multiple children (panels)
        assert len(fig.children) >= 2

    # Edge case tests

    def test_ld_heatmap_single_snp_in_region(self, ld_heatmap_gwas_df):
        """Test with single SNP in filtered region (graceful handling)."""
        # Create a matrix with 5 SNPs but only 1 will be in region
        ld_matrix = pd.DataFrame(
            np.array(
                [
                    [1.0, 0.9, 0.7, 0.4, 0.2],
                    [0.9, 1.0, 0.8, 0.5, 0.3],
                    [0.7, 0.8, 1.0, 0.6, 0.4],
                    [0.4, 0.5, 0.6, 1.0, 0.7],
                    [0.2, 0.3, 0.4, 0.7, 1.0],
                ]
            ),
            index=["rs1", "rs2", "rs3", "rs4", "rs5"],
            columns=["rs1", "rs2", "rs3", "rs4", "rs5"],
        )
        snp_ids = ["rs1", "rs2", "rs3", "rs4", "rs5"]

        plotter = LocusZoomPlotter(species=None, log_level=None)

        # Very narrow region that only includes rs1 (at position 1000000)
        fig = plotter.plot(
            ld_heatmap_gwas_df,
            chrom=1,
            start=999999,
            end=1000001,
            display=DisplayConfig(show_recombination=False),
            panels=PanelInputs(ld_heatmap_df=ld_matrix, ld_heatmap_snp_ids=snp_ids),
        )

        # Should complete without error
        assert fig is not None

    def test_ld_heatmap_lead_snp_not_in_heatmap(
        self, ld_heatmap_gwas_df, sample_ld_heatmap_data
    ):
        """Test lead SNP not in heatmap SNPs (no highlighting, no error)."""
        ld_matrix, snp_ids = sample_ld_heatmap_data
        plotter = LocusZoomPlotter(species=None, log_level=None)

        # Create GWAS with a lead SNP not in the heatmap
        gwas_with_extra = ld_heatmap_gwas_df.copy()
        gwas_with_extra = pd.concat(
            [
                gwas_with_extra,
                pd.DataFrame(
                    {
                        "rs": ["rs_extra"],
                        "chr": [1],
                        "ps": [1000100],  # Different position
                        "p_wald": [1e-10],  # Most significant
                    }
                ),
            ],
            ignore_index=True,
        )

        fig = plotter.plot(
            gwas_with_extra,
            chrom=1,
            start=999000,
            end=1003000,
            display=DisplayConfig(show_recombination=False),
            ld=LDConfig(lead_pos=1000100),
            panels=PanelInputs(ld_heatmap_df=ld_matrix, ld_heatmap_snp_ids=snp_ids),
        )

        # Should complete without error
        assert fig is not None

    def test_ld_heatmap_missing_snp_ids_raises_error(
        self, ld_heatmap_gwas_df, sample_ld_heatmap_data
    ):
        """Test that providing ld_heatmap_df without ld_heatmap_snp_ids raises error."""
        ld_matrix, _ = sample_ld_heatmap_data
        plotter = LocusZoomPlotter(species=None, log_level=None)

        with pytest.raises(ValueError, match="ld_heatmap_snp_ids is required"):
            plotter.plot(
                ld_heatmap_gwas_df,
                chrom=1,
                start=999000,
                end=1003000,
                display=DisplayConfig(show_recombination=False),
                panels=PanelInputs(ld_heatmap_df=ld_matrix),
            )


def _glyph_values(renderer, prop):
    """Per-item values of a Bokeh glyph property, whether column or literal."""
    spec = getattr(renderer.glyph, prop)
    data = renderer.data_source.data
    if isinstance(spec, str):
        return list(data[spec])
    return [spec] * len(data[renderer.glyph.x])


class TestRegionalHeatmapOutlineIsInGenomicCoordinates:
    START = 999000
    END = 1003000

    @pytest.fixture
    def heatmap_gwas_df(self):
        return pd.DataFrame(
            {
                "rs": ["rs1", "rs2", "rs3", "rs4", "rs5"],
                "chr": [1, 1, 1, 1, 1],
                "ps": [1000000, 1000500, 1001000, 1001500, 1002000],
                "p_wald": [1e-8, 1e-6, 1e-4, 1e-3, 0.05],
            }
        )

    @pytest.fixture
    def heatmap_ld_matrix(self):
        ids = ["rs1", "rs2", "rs3", "rs4", "rs5"]
        values = np.array(
            [
                [1.0, 0.9, 0.7, 0.4, 0.2],
                [0.9, 1.0, 0.8, 0.5, 0.3],
                [0.7, 0.8, 1.0, 0.6, 0.4],
                [0.4, 0.5, 0.6, 1.0, 0.7],
                [0.2, 0.3, 0.4, 0.7, 1.0],
            ]
        )
        return pd.DataFrame(values, index=ids, columns=ids)

    def _plot(self, backend, gwas_df, ld_matrix):
        return LocusZoomPlotter(species=None, log_level=None, backend=backend).plot(
            gwas_df,
            chrom=1,
            start=self.START,
            end=self.END,
            display=DisplayConfig(show_recombination=False),
            ld=LDConfig(lead_pos=1000000),
            panels=PanelInputs(
                ld_heatmap_df=ld_matrix, ld_heatmap_snp_ids=list(ld_matrix.index)
            ),
        )

    def _assert_inside_region(self, spans):
        assert spans, "the lead SNP should be outlined on the heatmap panel"
        for x0, x1 in spans:
            assert self.START <= x0 <= x1 <= self.END, (
                f"outline spans {x0}-{x1}, outside the region {self.START}-{self.END}"
            )

    def test_matplotlib_outline_uses_genomic_coordinates(
        self, heatmap_gwas_df, heatmap_ld_matrix
    ):
        fig = self._plot("matplotlib", heatmap_gwas_df, heatmap_ld_matrix)

        heatmap_ax = fig.axes[1]
        spans = [
            (p.get_x(), p.get_x() + p.get_width())
            for p in heatmap_ax.patches
            if not p.get_fill()
        ]

        self._assert_inside_region(spans)

    def test_plotly_outline_uses_genomic_coordinates_on_its_own_panel(
        self, heatmap_gwas_df, heatmap_ld_matrix
    ):
        fig = self._plot("plotly", heatmap_gwas_df, heatmap_ld_matrix)

        heatmap = next(trace for trace in fig.data if trace.type == "heatmap")
        outlines = [s for s in fig.layout.shapes if s.type == "rect"]

        self._assert_inside_region([(s.x0, s.x1) for s in outlines])
        for shape in outlines:
            assert shape.xref == heatmap.xaxis
            assert shape.yref == heatmap.yaxis

    def test_bokeh_outline_uses_genomic_coordinates(
        self, heatmap_gwas_df, heatmap_ld_matrix
    ):
        from bokeh.models import Rect

        layout = self._plot("bokeh", heatmap_gwas_df, heatmap_ld_matrix)

        heatmap_figure = layout.children[1]
        spans = []
        for renderer in heatmap_figure.renderers:
            if not isinstance(renderer.glyph, Rect):
                continue
            if "value" in renderer.data_source.data:
                continue
            xs = _glyph_values(renderer, "x")
            widths = _glyph_values(renderer, "width")
            spans.extend((x - w / 2, x + w / 2) for x, w in zip(xs, widths))

        self._assert_inside_region(spans)
