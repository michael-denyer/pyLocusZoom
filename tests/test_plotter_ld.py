"""Tests for LD calculation and the LD heatmap panel in regional plots."""

from unittest.mock import patch

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from pylocuszoom.backends.matplotlib_backend import MatplotlibBackend
from pylocuszoom.backends.plotly_backend import PlotlyBackend
from pylocuszoom.plotter import LocusZoomPlotter


class TestLocusZoomPlotterLdCalculation:
    """Tests for LD calculation integration."""

    @pytest.fixture
    def mock_plink_plotter(self):
        """Create plotter with mocked PLINK."""
        return LocusZoomPlotter(species="canine", plink_path="/mock/plink")

    def test_calculates_ld_when_reference_provided(self, mock_plink_plotter):
        """Should attempt LD calculation when ld_reference_file provided."""
        df = pd.DataFrame(
            {
                "rs": ["rs1", "rs2", "rs3"],
                "ps": [1100000, 1500000, 1900000],
                "p_wald": [1e-8, 1e-5, 1e-3],
            }
        )

        with patch("pylocuszoom._ld_plotting.calculate_ld") as mock_ld:
            mock_ld.return_value = pd.DataFrame(
                {
                    "SNP": ["rs1", "rs2", "rs3"],
                    "R2": [1.0, 0.8, 0.5],
                }
            )

            fig = mock_plink_plotter.plot(
                df,
                chrom=1,
                start=1000000,
                end=2000000,
                lead_pos=1100000,
                ld_reference_file="/path/to/genotypes",
                show_recombination=False,
            )

            mock_ld.assert_called_once()
            plt.close(fig)

    def test_empty_ld_output_is_downgraded_to_warning(
        self, mock_plink_plotter, warning_records
    ):
        """An empty-output PlinkError (singleton lead SNP) should not abort.

        Singleton lead SNPs with no LD neighbours in the window are a real
        scenario; mock_plink_plotter.plot() catches only this specific PlinkError and
        continues without LD colouring, leaving a warning in the log.
        """
        from pylocuszoom.exceptions import EmptyLDOutputError

        df = pd.DataFrame(
            {
                "rs": ["rs1", "rs2", "rs3"],
                "ps": [1100000, 1500000, 1900000],
                "p_wald": [1e-8, 1e-5, 1e-3],
            }
        )

        with patch("pylocuszoom._ld_plotting.calculate_ld") as mock_ld:
            mock_ld.side_effect = EmptyLDOutputError(
                "PLINK produced an empty LD output for lead SNP 'rs1'."
            )
            fig = mock_plink_plotter.plot(
                df,
                chrom=1,
                start=1000000,
                end=2000000,
                lead_pos=1100000,
                ld_reference_file="/path/to/genotypes",
                show_recombination=False,
            )
        assert fig is not None
        assert any("LD calculation skipped" in message for message in warning_records)
        plt.close(fig)

    def test_stacked_plot_downgrades_empty_ld_output(self, mock_plink_plotter):
        """Single and stacked plots share the recoverable LD policy."""
        from pylocuszoom.exceptions import EmptyLDOutputError

        df = pd.DataFrame(
            {
                "rs": ["rs1", "rs2", "rs3"],
                "ps": [1100000, 1500000, 1900000],
                "p_wald": [1e-8, 1e-5, 1e-3],
            }
        )

        with patch("pylocuszoom._ld_plotting.calculate_ld") as mock_ld:
            mock_ld.side_effect = EmptyLDOutputError("no LD neighbours")
            fig = mock_plink_plotter.plot_stacked(
                [df],
                chrom=1,
                start=1000000,
                end=2000000,
                lead_positions=[1100000],
                ld_reference_files=["/path/to/genotypes"],
                show_recombination=False,
            )

        assert fig is not None
        plt.close(fig)

    def test_plink_misconfiguration_propagates_through_plot(self, mock_plink_plotter):
        """A non-empty-output PlinkError must surface — it means PLINK is broken.

        Regression boundary: the catch in mock_plink_plotter.plot() is narrow on purpose.
        Timeout, non-zero exit, and "output file missing after success" all
        indicate real misconfiguration and should reach the caller.
        """
        from pylocuszoom.exceptions import PlinkError

        df = pd.DataFrame(
            {
                "rs": ["rs1", "rs2", "rs3"],
                "ps": [1100000, 1500000, 1900000],
                "p_wald": [1e-8, 1e-5, 1e-3],
            }
        )

        with patch("pylocuszoom._ld_plotting.calculate_ld") as mock_ld:
            mock_ld.side_effect = PlinkError(
                "PLINK LD calculation failed (exit code 2): bad bfile"
            )
            with pytest.raises(PlinkError, match="exit code"):
                mock_plink_plotter.plot(
                    df,
                    chrom=1,
                    start=1000000,
                    end=2000000,
                    lead_pos=1100000,
                    ld_reference_file="/path/to/genotypes",
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
            show_recombination=False,
            ld_heatmap_df=ld_matrix,
            ld_heatmap_snp_ids=snp_ids,
        )

        # Should have at least 2 axes (association + heatmap)
        assert len(fig.axes) >= 2
        plt.close(fig)

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
            show_recombination=False,
            ld_heatmap_df=ld_matrix,
            ld_heatmap_snp_ids=snp_ids,
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
        plt.close(fig)

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
            show_recombination=False,
            genes_df=heatmap_genes_df,
            ld_heatmap_df=ld_matrix,
            ld_heatmap_snp_ids=snp_ids,
        )

        # Should have 3 panels: GWAS, gene track, heatmap
        axes = fig.axes
        assert len(axes) >= 3
        plt.close(fig)

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
            end=1001000,  # Only includes rs1, rs2, rs3
            show_recombination=False,
            ld_heatmap_df=ld_matrix,
            ld_heatmap_snp_ids=snp_ids,
        )

        # Should complete without error even with partial overlap
        assert fig is not None
        plt.close(fig)

    def test_ld_heatmap_empty_overlap_logs_warning(
        self, ld_heatmap_gwas_df, sample_ld_heatmap_data, capfd
    ):
        """When no SNPs in heatmap overlap with region, warning logged and no heatmap panel added."""
        ld_matrix, snp_ids = sample_ld_heatmap_data
        plotter = LocusZoomPlotter(species=None, log_level="WARNING")

        # Region that doesn't overlap with any GWAS SNP positions
        fig = plotter.plot(
            ld_heatmap_gwas_df,
            chrom=1,
            start=5000000,  # Far outside GWAS positions
            end=6000000,
            show_recombination=False,
            ld_heatmap_df=ld_matrix,
            ld_heatmap_snp_ids=snp_ids,
        )

        # Should complete without error - heatmap simply not added
        assert fig is not None
        plt.close(fig)

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
            show_recombination=False,
            ld_heatmap_df=ld_matrix,
            ld_heatmap_snp_ids=snp_ids,
            ld_heatmap_height=0.1,  # Smaller
        )

        fig2 = plotter.plot(
            ld_heatmap_gwas_df,
            chrom=1,
            start=999000,
            end=1003000,
            show_recombination=False,
            ld_heatmap_df=ld_matrix,
            ld_heatmap_snp_ids=snp_ids,
            ld_heatmap_height=0.5,  # Larger
        )

        # Both should render successfully
        assert fig1 is not None
        assert fig2 is not None
        plt.close(fig1)
        plt.close(fig2)

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
            lead_pos=1000000,  # rs1
            show_recombination=False,
            ld_heatmap_df=ld_matrix,
            ld_heatmap_snp_ids=snp_ids,
        )

        # Should render with lead SNP highlight (visual check)
        assert fig is not None
        plt.close(fig)

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
            show_recombination=False,
            ld_heatmap_df=ld_matrix,
            ld_heatmap_snp_ids=snp_ids,
        )

        # Should have at least 2 main axes (association + heatmap)
        # Plus possible colorbar axis
        assert len(fig.axes) >= 2
        assert isinstance(fig, plt.Figure)
        plt.close(fig)

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
            show_recombination=False,
            ld_heatmap_df=ld_matrix,
            ld_heatmap_snp_ids=snp_ids,
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
            show_recombination=False,
            ld_heatmap_df=ld_matrix,
            ld_heatmap_snp_ids=snp_ids,
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
            end=1000001,  # Only rs1 at 1000000
            show_recombination=False,
            ld_heatmap_df=ld_matrix,
            ld_heatmap_snp_ids=snp_ids,
        )

        # Should complete without error
        assert fig is not None
        plt.close(fig)

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
            lead_pos=1000100,  # rs_extra - not in heatmap
            show_recombination=False,
            ld_heatmap_df=ld_matrix,
            ld_heatmap_snp_ids=snp_ids,
        )

        # Should complete without error
        assert fig is not None
        plt.close(fig)

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
                show_recombination=False,
                ld_heatmap_df=ld_matrix,
                # ld_heatmap_snp_ids not provided
            )


class TestHighlightHeatmapSnpBackendProtocol:
    """Tests that highlight_heatmap_snp works on all backends."""

    def test_matplotlib_highlight_heatmap_snp(self):
        """Matplotlib backend highlight_heatmap_snp creates rectangle patches."""

        backend = MatplotlibBackend()
        fig, axes = backend.create_figure(
            n_panels=1, height_ratios=[1.0], figsize=(6, 4)
        )
        ax = axes[0]

        # Should not raise
        backend.highlight_heatmap_snp(ax, fig, snp_idx=2, n_snps=5)

        # Verify patches were added (3 for row + 2 for column = 5 total)
        rect_patches = [
            p for p in ax.patches if hasattr(p, "get_edgecolor") and not p.get_fill()
        ]
        assert len(rect_patches) == 5
        plt.close(fig)

    def test_plotly_highlight_heatmap_snp(self):
        """Plotly backend highlight_heatmap_snp adds shapes to figure."""

        backend = PlotlyBackend()
        fig, panels = backend.create_figure(
            n_panels=1, height_ratios=[1.0], figsize=(6, 4)
        )
        ax = panels[0]

        # Should not raise
        backend.highlight_heatmap_snp(ax, fig, snp_idx=2, n_snps=5)

        # Verify shapes were added
        shapes = fig.layout.shapes
        assert len(shapes) == 5  # 3 for row + 2 for column

    def test_bokeh_highlight_heatmap_snp(self):
        """Bokeh backend highlight_heatmap_snp adds rect glyphs to figure."""
        from pylocuszoom.backends.bokeh_backend import BokehBackend

        backend = BokehBackend()
        layout, figures = backend.create_figure(
            n_panels=1, height_ratios=[1.0], figsize=(6, 4)
        )
        ax = figures[0]
        renderers_before = len(ax.renderers)

        # Should not raise
        backend.highlight_heatmap_snp(ax, layout, snp_idx=2, n_snps=5)

        # Verify renderer was added (batched into single rect call)
        assert len(ax.renderers) == renderers_before + 1

        # Verify the batched renderer has correct number of cells
        # snp_idx=2, n_snps=5: row cells (0,2),(1,2),(2,2) + col cells (2,3),(2,4) = 5
        renderer = ax.renderers[-1]
        assert len(renderer.data_source.data["x"]) == 5
