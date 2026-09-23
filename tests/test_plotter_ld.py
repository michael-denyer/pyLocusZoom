"""Tests for LD calculation and the LD heatmap panel in regional plots."""

import numpy as np
import pandas as pd
import pytest

from pylocuszoom import DisplayConfig, LDConfig, PanelInputs
from pylocuszoom.backends import BUILTIN_BACKENDS
from pylocuszoom.backends.composition import LD_LEGEND_TITLE
from pylocuszoom.exceptions import PlinkError
from pylocuszoom.plotter import LocusZoomPlotter
from tests.conftest import FIGURE_TYPES
from tests.figure_probes import PROBES


def _drawn_positions(ax):
    """Return the distinct x coordinates of every point scattered on ``ax``."""
    return {
        float(x) for collection in ax.collections for x, _ in collection.get_offsets()
    }


def _heatmap_height_ratio(fig):
    """Return the heatmap panel's drawn height as a fraction of the association panel's."""
    association, heatmap = fig.get_axes()[:2]
    return heatmap.get_position().height / association.get_position().height


class TestLocusZoomPlotterLdCalculation:
    """Tests for LD calculation integration."""

    LD_OUTPUT = (
        "CHR_A BP_A SNP_A CHR_B BP_B SNP_B R2\n"
        "1 1100000 rs1 1 1500000 rs2 0.80\n"
        "1 1100000 rs1 1 1900000 rs3 0.50\n"
    )

    @pytest.mark.parametrize("existing_r2", [False, True])
    def test_ld_reference_colours_the_plot(
        self, fake_plink, tiny_regional_gwas_df, existing_r2
    ):
        """A PLINK run that returns LD pairs colours the points by R2."""
        bfile, plink_writes = fake_plink
        frame = (
            tiny_regional_gwas_df.assign(R2=0.0)
            if existing_r2
            else tiny_regional_gwas_df
        )
        original = frame.copy(deep=True)

        with plink_writes(self.LD_OUTPUT):
            fig = LocusZoomPlotter(species="canine", plink_path="/mock/plink").plot(
                frame,
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
        from matplotlib.colors import to_hex

        colors = {
            to_hex(points.get_facecolors()[0]) for points in fig.axes[0].collections
        }
        assert {"#ff0000", "#00cd00"} <= colors
        pd.testing.assert_frame_equal(frame, original)

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
        """Draw the stacked panel with no r² legend when PLINK returns no pairs."""
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

        (ax,) = fig.get_axes()
        assert ax.get_legend() is None, "no LD pairs means no r² legend"
        assert _drawn_positions(ax) == {1100000.0, 1500000.0, 1900000.0}

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
                display=DisplayConfig(show_recombination=False),
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
                "pos": [1000000, 1000500, 1001000, 1001500, 1002000],
                "p_value": [1e-8, 1e-6, 1e-4, 1e-3, 0.05],
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
        """Keep only the heatmap SNPs inside [start, end], shrinking the matrix."""
        ld_matrix, snp_ids = sample_ld_heatmap_data
        plotter = LocusZoomPlotter(species=None, log_level=None)

        fig = plotter.plot(
            ld_heatmap_gwas_df,
            chrom=1,
            start=1000000,
            end=1001000,
            display=DisplayConfig(show_recombination=False),
            panels=PanelInputs(ld_heatmap_df=ld_matrix, ld_heatmap_snp_ids=snp_ids),
        )

        (image,) = fig.get_axes()[1].collections
        assert image.get_array().shape == (3, 3), "rs1, rs2 and rs3 are in the region"
        assert image.get_coordinates()[0, :, 0].tolist() == [
            999750,
            1000250,
            1000750,
            1001250,
        ]

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
        """Size the heatmap panel as ld_heatmap_height of the association panel."""
        ld_matrix, snp_ids = sample_ld_heatmap_data
        plotter = LocusZoomPlotter(species=None, log_level=None)

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

        assert _heatmap_height_ratio(fig1) == pytest.approx(0.1, abs=1e-3)
        assert _heatmap_height_ratio(fig2) == pytest.approx(0.5, abs=1e-3)

    def test_ld_heatmap_lead_snp_highlight(
        self, ld_heatmap_gwas_df, sample_ld_heatmap_data
    ):
        """Draw the highlight crosshair over the lead SNP's own heatmap column."""
        ld_matrix, snp_ids = sample_ld_heatmap_data
        plotter = LocusZoomPlotter(species=None, log_level=None)

        fig = plotter.plot(
            ld_heatmap_gwas_df,
            chrom=1,
            start=999000,
            end=1003000,
            display=DisplayConfig(show_recombination=False),
            ld=LDConfig(lead_pos=1000000),
            panels=PanelInputs(ld_heatmap_df=ld_matrix, ld_heatmap_snp_ids=snp_ids),
        )

        highlights = fig.get_axes()[1].patches
        assert highlights, "rs1 is in the heatmap, so it gets a highlight"
        assert {(rect.get_x(), rect.get_width()) for rect in highlights} == {
            (999750.0, 500.0)
        }

    @pytest.mark.parametrize("backend_name", BUILTIN_BACKENDS)
    def test_ld_heatmap_adds_a_panel_under_the_association_plot(
        self, backend_name, ld_heatmap_gwas_df, sample_ld_heatmap_data
    ):
        """Every backend stacks the heatmap as a second panel."""
        ld_matrix, snp_ids = sample_ld_heatmap_data
        plotter = LocusZoomPlotter(species=None, backend=backend_name, log_level=None)

        fig = plotter.plot(
            ld_heatmap_gwas_df,
            chrom=1,
            start=999000,
            end=1003000,
            display=DisplayConfig(show_recombination=False),
            panels=PanelInputs(ld_heatmap_df=ld_matrix, ld_heatmap_snp_ids=snp_ids),
        )

        assert isinstance(fig, FIGURE_TYPES[backend_name])
        assert PROBES[backend_name].panel_count(fig) == 2

    # Edge case tests

    def test_ld_heatmap_single_snp_in_region(self, ld_heatmap_gwas_df):
        """Skip the heatmap, and its colorbar, when one SNP survives filtering."""
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

        fig = plotter.plot(
            ld_heatmap_gwas_df,
            chrom=1,
            start=999999,
            end=1000001,
            display=DisplayConfig(show_recombination=False),
            panels=PanelInputs(ld_heatmap_df=ld_matrix, ld_heatmap_snp_ids=snp_ids),
        )

        association, heatmap = fig.get_axes()
        assert association.collections, "the association panel still renders"
        assert len(heatmap.collections) == 0

    def test_ld_heatmap_lead_snp_not_in_heatmap(
        self, ld_heatmap_gwas_df, sample_ld_heatmap_data
    ):
        """Draw the heatmap unhighlighted when the lead SNP is not one of its SNPs."""
        ld_matrix, snp_ids = sample_ld_heatmap_data
        plotter = LocusZoomPlotter(species=None, log_level=None)

        gwas_with_extra = ld_heatmap_gwas_df.copy()
        gwas_with_extra = pd.concat(
            [
                gwas_with_extra,
                pd.DataFrame(
                    {
                        "rs": ["rs_extra"],
                        "chr": [1],
                        "pos": [1000100],  # Different position
                        "p_value": [1e-10],  # Most significant
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

        heatmap = fig.get_axes()[1]
        assert heatmap.collections, "the heatmap still renders"
        assert len(heatmap.patches) == 0, "rs_extra has no column to highlight"

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


class TestRegionalHeatmapOutlineIsInGenomicCoordinates:
    START = 999000
    END = 1003000

    @pytest.fixture
    def heatmap_gwas_df(self):
        return pd.DataFrame(
            {
                "rs": ["rs1", "rs2", "rs3", "rs4", "rs5"],
                "chr": [1, 1, 1, 1, 1],
                "pos": [1000000, 1000500, 1001000, 1001500, 1002000],
                "p_value": [1e-8, 1e-6, 1e-4, 1e-3, 0.05],
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

    @pytest.mark.parametrize("backend_name", BUILTIN_BACKENDS)
    def test_outline_uses_genomic_coordinates_on_the_heatmap_panel(
        self, backend_name, heatmap_gwas_df, heatmap_ld_matrix
    ):
        fig = self._plot(backend_name, heatmap_gwas_df, heatmap_ld_matrix)

        outlines = [
            b for b in PROBES[backend_name].boxes(fig, 1) if b.facecolor is None
        ]

        self._assert_inside_region([(b.x0, b.x1) for b in outlines])


def test_regional_heatmap_sorts_coordinates_and_matrix_together():
    frame = pd.DataFrame(
        {
            "chr": 1,
            "pos": [100, 200, 1000],
            "p_value": [0.1, 0.01, 0.001],
            "rs": ["a", "b", "c"],
        }
    )
    matrix = pd.DataFrame([[0.9, 0.1, 0.2], [0.1, 0.8, 0.3], [0.2, 0.3, 0.7]])
    fig = LocusZoomPlotter(species=None, backend="bokeh", log_level=None).plot(
        frame,
        chrom=1,
        start=1,
        end=1500,
        display=DisplayConfig(show_recombination=False, snp_labels=False),
        panels=PanelInputs(ld_heatmap_df=matrix, ld_heatmap_snp_ids=["c", "a", "b"]),
    )
    cells = fig.children[-1].renderers[0].data_source.data
    bounds = [(x - w / 2, x + w / 2) for x, w in zip(cells["x"], cells["w"])]
    assert bounds == [
        (50, 150),
        (50, 150),
        (150, 600),
        (50, 150),
        (150, 600),
        (600, 1400),
    ]
    assert cells["value"] == [0.8, 0.3, 0.7, 0.1, 0.2, 0.9]


def test_regional_heatmap_rejects_duplicate_genomic_coordinates():
    frame = pd.DataFrame(
        {"chr": 1, "pos": [150, 150], "p_value": [0.1, 0.01], "rs": ["a", "b"]}
    )
    with pytest.raises(ValueError, match="distinct genomic positions"):
        LocusZoomPlotter(species=None, log_level=None).plot(
            frame,
            chrom=1,
            start=100,
            end=200,
            display=DisplayConfig(show_recombination=False, snp_labels=False),
            panels=PanelInputs(
                ld_heatmap_df=pd.DataFrame(np.eye(2)), ld_heatmap_snp_ids=["a", "b"]
            ),
        )


def test_heatmap_highlights_selected_variant_at_duplicate_source_position():
    frame = pd.DataFrame(
        {
            "chr": 1,
            "pos": [150, 150, 250],
            "p_value": [0.1, 1e-8, 0.001],
            "rs": ["weak", "strong", "other"],
        }
    )
    fig = LocusZoomPlotter(species=None, log_level=None).plot(
        frame,
        chrom=1,
        start=100,
        end=300,
        display=DisplayConfig(show_recombination=False, snp_labels=False),
        panels=PanelInputs(
            ld_heatmap_df=pd.DataFrame(np.eye(2)),
            ld_heatmap_snp_ids=["other", "strong"],
        ),
    )
    outlines = fig.axes[1].patches
    assert [patch.get_xy() for patch in outlines] == [(100, -0.5), (100, 0.5)]
    assert [patch.get_width() for patch in outlines] == [100, 100]


@pytest.mark.parametrize("with_genes", [False, True])
@pytest.mark.parametrize("with_recombination", [False, True])
def test_regional_colorbar_preserves_genomic_display_alignment(
    sample_genes_df, sample_recomb_df, with_genes, with_recombination
):
    frame = pd.DataFrame(
        {
            "chr": 1,
            "pos": [1_100_000, 1_200_000, 1_900_000],
            "p_value": [0.1, 0.01, 0.001],
            "rs": ["a", "b", "c"],
        }
    )
    fig = LocusZoomPlotter(species=None, log_level=None).plot(
        frame,
        chrom=1,
        start=1_000_000,
        end=2_000_000,
        display=DisplayConfig(show_recombination=with_recombination, snp_labels=False),
        panels=PanelInputs(
            ld_heatmap_df=pd.DataFrame(np.eye(3)),
            ld_heatmap_snp_ids=["a", "b", "c"],
            genes_df=sample_genes_df if with_genes else None,
            recomb_df=sample_recomb_df if with_recombination else None,
        ),
    )
    fig.canvas.draw()
    association = fig.axes[0]
    for panel in association.get_shared_x_axes().get_siblings(association):
        np.testing.assert_allclose(
            panel.transData.transform([[1_100_000, 0], [1_900_000, 0]])[:, 0],
            association.transData.transform([[1_100_000, 0], [1_900_000, 0]])[:, 0],
        )
    colorbar = next(axis for axis in fig.axes if axis.get_ylabel() == "R²")
    assert colorbar.get_position().x0 > association.get_position().x1
    assert colorbar.get_tightbbox(fig.canvas.get_renderer()).x1 <= fig.bbox.x1


class TestPlotEdgeCases:
    """Tests for plot() edge cases and error handling."""

    @pytest.fixture
    def mock_plink_plotter(self):
        """Create plotter instance."""
        return LocusZoomPlotter(species="canine", plink_path="/mock/plink")

    def test_plot_skips_ld_when_rs_col_missing(
        self, mock_plink_plotter, warning_records
    ):
        """LD needs SNP IDs, so a frame without rs_col plots uncoloured.

        The alternative, a KeyError from deep inside the LD merge, tells the
        caller nothing about which column is missing.
        """
        df = pd.DataFrame(
            {
                "chr": 1,
                "pos": [1100000, 1500000, 1900000],
                "p_value": [1e-8, 1e-5, 1e-3],
            }
        )

        fig = mock_plink_plotter.plot(
            df,
            chrom=1,
            start=1000000,
            end=2000000,
            display=DisplayConfig(show_recombination=False),
            ld=LDConfig(lead_pos=1500000, ld_reference_file="/path/to/genotypes"),
        )

        assert fig.get_axes()[0].get_legend() is None
        assert any("rs" in message for message in warning_records)
