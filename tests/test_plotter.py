"""Tests for LocusZoomPlotter class."""

from pathlib import Path
from unittest.mock import Mock, patch

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from matplotlib.figure import Figure

from pylocuszoom._gene_source import GeneAnnotations
from pylocuszoom._regional_panels import AssociationPanel
from pylocuszoom.plotter import LocusZoomPlotter
from tests.reference_mocks import (
    gene_transcript_exon_payload,
    ok_response,
    refseq_payload,
)


class TestLocusZoomPlotterInit:
    """Tests for LocusZoomPlotter initialization."""

    def test_default_species_is_canine(self):
        """Default species should be canine."""
        plotter = LocusZoomPlotter()
        assert plotter.species == "canine"

    def test_custom_species(self):
        """Should accept custom species."""
        plotter = LocusZoomPlotter(species="feline")
        assert plotter.species == "feline"

    def test_custom_plink_path(self):
        """Should accept custom PLINK path."""
        plotter = LocusZoomPlotter(plink_path="/custom/plink")
        assert plotter.plink_path == "/custom/plink"

    def test_custom_threshold(self):
        """Should accept custom genomewide threshold."""
        plotter = LocusZoomPlotter(genomewide_threshold=5e-8)
        assert plotter.genomewide_threshold == 5e-8

    def test_auto_genes_default_false(self):
        """auto_genes should be False by default for backward compatibility."""
        plotter = LocusZoomPlotter(species="canine", log_level=None)
        assert plotter._auto_genes is False

    def test_auto_genes_can_be_enabled(self):
        """auto_genes=True should be accepted."""
        plotter = LocusZoomPlotter(species="human", log_level=None, auto_genes=True)
        assert plotter._auto_genes is True


class TestAutoGenes:
    """Tests for automatic gene fetching from Ensembl."""

    @pytest.fixture
    def two_gene_track_df(self):
        """Sample gene DataFrame for testing."""
        return pd.DataFrame(
            {
                "chr": ["1", "1"],
                "start": [1100000, 1300000],
                "end": [1200000, 1400000],
                "gene_name": ["GENE1", "GENE2"],
                "strand": ["+", "-"],
            }
        )

    def test_plot_with_auto_genes_enabled(self, small_regional_gwas_df):
        """Test that auto_genes=True fetches genes from the reference source."""
        # Mock the Ensembl API response
        mock_genes = pd.DataFrame(
            {
                "chr": ["1", "1"],
                "start": [1000000, 1500000],
                "end": [1200000, 1700000],
                "gene_name": ["GENE1", "GENE2"],
                "strand": ["+", "-"],
            }
        )

        plotter = LocusZoomPlotter(species="human", log_level=None, auto_genes=True)

        with patch(
            "pylocuszoom.plotter.get_genes_for_build",
            return_value=GeneAnnotations(mock_genes, pd.DataFrame()),
        ):
            fig = plotter.plot(
                small_regional_gwas_df,
                chrom=1,
                start=1000000,
                end=2000000,
            )

        assert fig is not None

    @pytest.mark.parametrize(
        "species,payload",
        [
            ("human", gene_transcript_exon_payload()),
            ("canine", refseq_payload()),
        ],
    )
    def test_auto_genes_draws_exons(
        self, small_regional_gwas_df, tmp_path, species, payload
    ):
        """auto_genes draws exon structure, whichever source served the genes."""
        from matplotlib.patches import Rectangle

        from pylocuszoom.gene_track import INTRON_HEIGHT

        plotter = LocusZoomPlotter(species=species, log_level=None, auto_genes=True)

        with (
            patch("pylocuszoom.reference_genes.cache_root", return_value=tmp_path),
            patch("pylocuszoom._http.requests.get", return_value=ok_response(payload)),
        ):
            fig = plotter.plot(
                small_regional_gwas_df,
                chrom=1,
                start=1_000_000,
                end=2_000_000,
                show_recombination=False,
            )

        introns = [
            rect
            for ax in fig.get_axes()
            for rect in ax.patches
            if isinstance(rect, Rectangle)
            and abs(rect.get_height() - INTRON_HEIGHT) < 1e-9
        ]
        assert introns, "no intron line drawn, so the gene track has no exons"

    def test_plot_auto_genes_warns_when_source_fails(
        self, small_regional_gwas_df, tmp_path
    ):
        """A gene-source outage warns instead of passing as an empty region."""
        plotter = LocusZoomPlotter(species="canine", log_level=None, auto_genes=True)
        unavailable = Mock(ok=False, status_code=503, text="Service Unavailable")

        with (
            patch("pylocuszoom.reference_genes.cache_root", return_value=tmp_path),
            patch("pylocuszoom._http.time.sleep"),
            patch("pylocuszoom._http.requests.get", return_value=unavailable),
            pytest.warns(UserWarning, match=r"chr1:1000000-2000000.*UCSC.*503"),
        ):
            fig = plotter.plot(
                small_regional_gwas_df,
                chrom=1,
                start=1000000,
                end=2000000,
                show_recombination=False,
            )

        assert fig is not None

    def test_plot_stacked_auto_genes_overrides_constructor(
        self, small_regional_gwas_df
    ):
        """plot_stacked(auto_genes=True) fetches genes on a plotter built without it."""
        mock_genes = pd.DataFrame(
            {
                "chr": ["1"],
                "start": [1200000],
                "end": [1400000],
                "gene_name": ["GENE1"],
                "strand": ["+"],
            }
        )
        plotter = LocusZoomPlotter(species="human", log_level=None)

        with patch(
            "pylocuszoom.plotter.get_genes_for_build",
            return_value=GeneAnnotations(mock_genes, pd.DataFrame()),
        ) as mock_fetch:
            fig = plotter.plot_stacked(
                [small_regional_gwas_df],
                chrom=1,
                start=1000000,
                end=2000000,
                show_recombination=False,
                auto_genes=True,
            )

        mock_fetch.assert_called_once()
        assert len(fig.get_axes()) == 2

    def test_plot_auto_genes_disabled_by_default(self, small_regional_gwas_df):
        """Test that auto_genes=False by default (backward compatible)."""
        plotter = LocusZoomPlotter(species="canine", log_level=None)

        # Should work without genes_df and without calling Ensembl
        fig = plotter.plot(
            small_regional_gwas_df,
            chrom=1,
            start=1000000,
            end=2000000,
        )

        assert fig is not None

    def test_plot_auto_genes_respects_explicit_genes_df(
        self, small_regional_gwas_df, two_gene_track_df
    ):
        """Test that explicit genes_df is used even when auto_genes=True."""
        plotter = LocusZoomPlotter(species="human", log_level=None, auto_genes=True)

        with patch("pylocuszoom.plotter.get_genes_for_build") as mock_fetch:
            fig = plotter.plot(
                small_regional_gwas_df,
                chrom=1,
                start=1000000,
                end=2000000,
                genes_df=two_gene_track_df,
            )

            # Ensembl should NOT be called when genes_df is provided
            mock_fetch.assert_not_called()

        assert fig is not None


class TestLocusZoomPlotterPlot:
    """Tests for LocusZoomPlotter.plot() method."""

    def test_creates_figure(self, canine_plotter, regional_gwas_df):
        """Should create a matplotlib figure."""
        fig = canine_plotter.plot(
            regional_gwas_df,
            chrom=1,
            start=1000000,
            end=2000000,
        )
        assert isinstance(fig, plt.Figure)

    def test_plots_with_gene_track(
        self, canine_plotter, regional_gwas_df, sample_genes_df
    ):
        """Should create plot with gene track when genes_df provided."""
        fig = canine_plotter.plot(
            regional_gwas_df,
            chrom=1,
            start=1000000,
            end=2000000,
            genes_df=sample_genes_df,
        )
        # Should have 2 axes (association + gene track)
        assert len(fig.axes) >= 2

    def test_highlights_lead_snp(self, canine_plotter, regional_gwas_df):
        """Should highlight lead SNP when lead_pos provided."""
        lead_pos = regional_gwas_df["ps"].iloc[0]
        fig = canine_plotter.plot(
            regional_gwas_df,
            chrom=1,
            start=1000000,
            end=2000000,
            lead_pos=lead_pos,
        )
        assert isinstance(fig, Figure)
        # Test passes if no exception

    def test_handles_empty_dataframe(self, canine_plotter):
        """Empty GWAS DataFrame should raise ValidationError."""
        from pylocuszoom.exceptions import ValidationError

        empty_df = pd.DataFrame(columns=["rs", "chr", "ps", "p_wald"])
        with pytest.raises(ValidationError, match="empty"):
            canine_plotter.plot(
                empty_df,
                chrom=1,
                start=1000000,
                end=2000000,
            )

    def test_custom_column_names(self, canine_plotter):
        """Should work with custom column names."""
        df = pd.DataFrame(
            {
                "snp_id": ["rs1", "rs2", "rs3"],
                "position": [1100000, 1500000, 1900000],
                "pvalue": [1e-8, 1e-5, 1e-3],
            }
        )
        fig = canine_plotter.plot(
            df,
            chrom=1,
            start=1000000,
            end=2000000,
            pos_col="position",
            p_col="pvalue",
            rs_col="snp_id",
        )
        assert isinstance(fig, Figure)

    def test_with_precomputed_ld(self, canine_plotter, regional_gwas_df):
        """Should use pre-computed LD column when provided."""
        df = regional_gwas_df.copy()
        df["R2"] = np.random.default_rng(0).uniform(0, 1, len(df))

        fig = canine_plotter.plot(
            df,
            chrom=1,
            start=1000000,
            end=2000000,
            ld_col="R2",
        )
        assert isinstance(fig, Figure)

    def test_with_recombination_data(
        self, canine_plotter, regional_gwas_df, sample_recomb_df
    ):
        """Should plot with recombination overlay when provided."""
        fig = canine_plotter.plot(
            regional_gwas_df,
            chrom=1,
            start=1000000,
            end=2000000,
            recomb_df=sample_recomb_df,
        )
        assert isinstance(fig, Figure)

    def test_disables_snp_labels(self, canine_plotter, regional_gwas_df):
        """Should not add labels when snp_labels=False."""
        fig = canine_plotter.plot(
            regional_gwas_df,
            chrom=1,
            start=1000000,
            end=2000000,
            snp_labels=False,
        )
        assert isinstance(fig, Figure)

    def test_disables_recombination(self, canine_plotter, regional_gwas_df):
        """Should not show recombination when show_recombination=False."""
        fig = canine_plotter.plot(
            regional_gwas_df,
            chrom=1,
            start=1000000,
            end=2000000,
            show_recombination=False,
        )
        assert isinstance(fig, Figure)


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
                "ps": [1100000, 1500000, 1900000],
                "p_wald": [1e-8, 1e-5, 1e-3],
            }
        )

        fig = mock_plink_plotter.plot(
            df,
            chrom=1,
            start=1000000,
            end=2000000,
            lead_pos=1500000,
            ld_reference_file="/path/to/genotypes",
        )

        assert fig.get_axes()[0].get_legend() is None
        assert any("rs" in message for message in warning_records)


class TestPlotStackedEdgeCases:
    """Tests for plot_stacked() edge cases and error handling."""

    def test_plot_stacked_validates_eqtl_columns(
        self, canine_plotter, tiny_regional_gwas_df
    ):
        """Bug: plot_stacked() should validate eQTL DataFrame has required columns.

        Currently bypasses validate_eqtl_df() and directly accesses 'pos' and
        'p_value' columns at lines 945-952, causing cryptic KeyError instead
        of helpful validation message.
        """
        from pylocuszoom.eqtl import EQTLValidationError

        # eQTL data with wrong column names
        bad_eqtl_df = pd.DataFrame(
            {
                "position": [1500000],  # Should be 'pos'
                "pval": [1e-6],  # Should be 'p_value'
            }
        )

        # Should raise EQTLValidationError with helpful message
        # Currently raises KeyError: 'pos'
        with pytest.raises(EQTLValidationError):
            canine_plotter.plot_stacked(
                [tiny_regional_gwas_df],
                chrom=1,
                start=1000000,
                end=2000000,
                show_recombination=False,
                eqtl_df=bad_eqtl_df,
            )

    def test_plot_stacked_validates_list_lengths(
        self, canine_plotter, tiny_regional_gwas_df
    ):
        """Bug: plot_stacked() should error when list lengths don't match.

        Currently uses zip() which silently truncates the longer list.
        If user provides 3 GWAS DataFrames but only 2 lead_positions,
        the third GWAS is plotted without a lead SNP - confusing behavior.
        """
        gwas_dfs = [
            tiny_regional_gwas_df,
            tiny_regional_gwas_df.copy(),
            tiny_regional_gwas_df.copy(),
        ]
        lead_positions = [1500000, 1500000]  # Only 2, but 3 gwas_dfs

        # Should raise ValueError about mismatched lengths
        # Currently silently truncates - third GWAS has no lead SNP
        with pytest.raises(ValueError, match="lead_positions"):
            canine_plotter.plot_stacked(
                gwas_dfs,
                chrom=1,
                start=1000000,
                end=2000000,
                lead_positions=lead_positions,
                show_recombination=False,
            )

    def test_plot_stacked_validates_panel_labels_length(
        self, canine_plotter, tiny_regional_gwas_df
    ):
        """Bug: panel_labels length should match gwas_dfs length."""
        gwas_dfs = [tiny_regional_gwas_df, tiny_regional_gwas_df.copy()]
        panel_labels = ["Only One"]  # Should have 2 labels

        # Should raise ValueError about mismatched lengths
        # Currently silently ignores - second panel has no label
        with pytest.raises(ValueError, match="panel_labels"):
            canine_plotter.plot_stacked(
                gwas_dfs,
                chrom=1,
                start=1000000,
                end=2000000,
                panel_labels=panel_labels,
                show_recombination=False,
            )


class TestPlotterDelegation:
    """Tests for plotter delegation to specialized classes."""

    def test_finemapping_panel_renders_the_supplied_pips(self):
        """plot_stacked() draws the fine-mapping frame onto its own PIP panel.

        Asserts on the plotted points rather than on a call to the renderer, so
        a change to the internal dispatch path cannot break this test without
        changing what the reader sees.
        """
        plotter = LocusZoomPlotter(
            species="canine", backend="matplotlib", log_level=None
        )

        gwas_df = pd.DataFrame({"ps": [1000, 2000], "p_wald": [0.01, 0.001]})
        fm_df = pd.DataFrame({"pos": [1000, 2000], "pip": [0.5, 0.3], "cs": [1, 1]})

        with patch.object(plotter, "_get_recomb_for_region", return_value=None):
            fig = plotter.plot_stacked(
                [gwas_df],
                chrom=1,
                start=1,
                end=3000,
                finemapping_df=fm_df,
            )

        pip_axes = [ax for ax in fig.get_axes() if ax.get_ylabel() == "PIP"]
        assert len(pip_axes) == 1, "fine-mapping data should add one PIP panel"
        plotted = pip_axes[0].collections[0].get_offsets().tolist()
        assert plotted == [[1000.0, 0.5], [2000.0, 0.3]]


def test_plotter_uses_ensure_recomb_maps():
    """Test that LocusZoomPlotter._ensure_recomb_maps delegates to module function."""
    with patch(
        "pylocuszoom.plotter.ensure_recomb_maps",
        return_value=Path("/mock/recomb"),
    ) as mock_ensure:
        plotter = LocusZoomPlotter(species="canine", log_level=None)

        # Call the method that should delegate to ensure_recomb_maps
        result = plotter._ensure_recomb_maps()

        mock_ensure.assert_called_once_with(species="canine", data_dir=None)
        assert result == Path("/mock/recomb")


# =============================================================================
# Property-Based Tests (Hypothesis)
# =============================================================================


# =============================================================================
# LD Heatmap Integration Tests
# =============================================================================


class TestStackedPlotLeadDetectionCrossChrom:
    """Regression: plot_stacked() lead auto-detection must filter by chromosome.

    Pre-fix bug: when lead_positions is None, the auto-detect loop did
    `df[(df[pos_col] >= start) & (df[pos_col] <= end)]` with no chrom
    filter. On a multi-chromosome GWAS DataFrame, the strongest p-value
    in [start, end] could come from a different chromosome than the
    plotted region, anchoring the diamond marker to the wrong locus.
    """

    def test_lead_autodetect_filters_by_chrom(self, canine_plotter):
        """Lead position must come from the requested chromosome."""
        # Two chromosomes share a position range. The strongest p-value
        # is on chr2 at position 1_500_000, but we are plotting chr1.
        # Pre-fix: lead = chr2's strongest hit. Post-fix: lead = chr1's.
        gwas_df = pd.DataFrame(
            {
                "rs": ["rs1", "rs2", "rs3", "rs4"],
                "chr": [1, 1, 2, 2],
                "ps": [1_200_000, 1_800_000, 1_500_000, 1_900_000],
                "p_wald": [1e-5, 1e-3, 1e-12, 1e-10],
            }
        )

        # Patch the lead-aware code to capture what got computed.
        captured = {}
        composer = canine_plotter._regional_composer
        original = composer.render_panel

        def spy(panel, *args, **kwargs):
            if isinstance(panel, AssociationPanel):
                captured.setdefault("lead_positions", []).append(panel.lead_pos)
            return original(panel, *args, **kwargs)

        composer.render_panel = spy
        try:
            canine_plotter.plot_stacked(
                [gwas_df],
                chrom=1,
                start=1_000_000,
                end=2_000_000,
                pos_col="ps",
                p_col="p_wald",
                show_recombination=False,
            )
        finally:
            composer.render_panel = original

        assert captured["lead_positions"] == [1_200_000], (
            "Lead must be chr1's strongest hit (1_200_000), not chr2's (1_500_000)"
        )


class TestLeadPosBoundary:
    """Pins the lead_pos=1 boundary.

    Smallest valid position (``ge=1`` in config) reaches the association
    panel intact, and the public API rejects ``0`` (1-based genomic coords).
    """

    def test_lead_pos_one_reaches_plot_association(self):
        """lead_pos=1 (smallest valid position) reaches the association scatter."""
        plotter = LocusZoomPlotter(species="canine", log_level=None)
        gwas_df = pd.DataFrame(
            {
                "rs": ["rs_lead", "rs2", "rs3"],
                "ps": [1, 100_000, 200_000],
                "p_wald": [1e-8, 1e-5, 1e-3],
            }
        )

        captured = {}
        composer = plotter._regional_composer
        original = composer.render_panel

        def spy(panel, *args, **kwargs):
            if isinstance(panel, AssociationPanel):
                captured["lead_pos"] = panel.lead_pos
            return original(panel, *args, **kwargs)

        composer.render_panel = spy
        try:
            plotter.plot(
                gwas_df,
                chrom=1,
                start=1,
                end=300_000,
                pos_col="ps",
                p_col="p_wald",
                lead_pos=1,
                show_recombination=False,
            )
        finally:
            composer.render_panel = original

        assert captured["lead_pos"] == 1, (
            "lead_pos=1 must pass through; falsy-check regression would drop it to None"
        )

    def test_lead_pos_zero_rejected_at_api(self):
        """Public API enforces genomic coords are 1-based; lead_pos=0 rejected."""
        plotter = LocusZoomPlotter(species="canine", log_level=None)
        gwas_df = pd.DataFrame(
            {
                "rs": ["rs1", "rs2"],
                "ps": [100_000, 200_000],
                "p_wald": [1e-8, 1e-5],
            }
        )

        from pydantic import ValidationError as PydanticValidationError

        with pytest.raises(PydanticValidationError, match="greater than or equal to 1"):
            plotter.plot(
                gwas_df,
                chrom=1,
                start=1,
                end=300_000,
                pos_col="ps",
                p_col="p_wald",
                lead_pos=0,
                show_recombination=False,
            )


class TestEmptyDataFrames:
    """Test handling of empty DataFrames across plot types.

    Behavior varies by plot type:
    - Regional plots: Allow empty DataFrames, render empty axes
    - Manhattan plots: Raise ValueError (can't compute axis limits)
    - QQ plots: Raise ValueError (require valid p-values)
    """

    @pytest.fixture
    def empty_gwas_df(self):
        """Empty DataFrame with correct columns."""
        return pd.DataFrame(columns=["chrom", "ps", "p_wald", "rs"])

    def test_plot_with_empty_df_raises(self, regional_plotter, empty_gwas_df):
        """Regional plot with empty DataFrame raises ValidationError."""
        from pylocuszoom.exceptions import ValidationError

        with pytest.raises(ValidationError, match="empty"):
            regional_plotter.plot(
                empty_gwas_df,
                chrom=1,
                start=1000000,
                end=2000000,
                pos_col="ps",
                p_col="p_wald",
                show_recombination=False,
            )


class TestNaNPvalues:
    """Test handling of NaN p-values in plots."""

    def test_plot_with_some_nan_pvalues_succeeds(self, regional_plotter):
        """Regional plot with partial NaN p-values should work.

        Rows with NaN p-values are excluded from plotting but do not
        cause errors. A warning is logged (visible in stderr).
        """
        gwas_df = pd.DataFrame(
            {
                "rs": ["rs1", "rs2", "rs3", "rs4", "rs5"],
                "ps": [1100000, 1300000, 1500000, 1700000, 1900000],
                "p_wald": [1e-8, np.nan, 1e-5, np.nan, 0.01],
            }
        )

        fig = regional_plotter.plot(
            gwas_df,
            chrom=1,
            start=1000000,
            end=2000000,
            show_recombination=False,
        )
        assert fig is not None

    def test_plot_with_all_nan_pvalues_succeeds(self, regional_plotter):
        """Regional plot with all NaN p-values should render empty.

        The plot renders without error but contains no data points.
        This tests graceful degradation rather than raising an error.
        """
        gwas_df = pd.DataFrame(
            {
                "rs": ["rs1", "rs2", "rs3"],
                "ps": [1100000, 1500000, 1900000],
                "p_wald": [np.nan, np.nan, np.nan],
            }
        )

        fig = regional_plotter.plot(
            gwas_df,
            chrom=1,
            start=1000000,
            end=2000000,
            show_recombination=False,
        )
        assert fig is not None


class TestStackedPlotMismatchedLengths:
    """Test mismatched list lengths in stacked plots raise ValidationError."""

    def test_stacked_mismatched_lead_positions_raises(
        self, regional_plotter, tiny_regional_gwas_df
    ):
        """Mismatched lead_positions length raises ValueError."""
        gwas_dfs = [
            tiny_regional_gwas_df,
            tiny_regional_gwas_df.copy(),
            tiny_regional_gwas_df.copy(),
        ]
        lead_positions = [1500000, 1500000]  # Only 2, but 3 gwas_dfs

        with pytest.raises(ValueError, match="lead_positions"):
            regional_plotter.plot_stacked(
                gwas_dfs,
                chrom=1,
                start=1000000,
                end=2000000,
                lead_positions=lead_positions,
                show_recombination=False,
            )

    def test_stacked_mismatched_panel_labels_raises(
        self, regional_plotter, tiny_regional_gwas_df
    ):
        """Mismatched panel_labels length raises ValueError."""
        gwas_dfs = [tiny_regional_gwas_df, tiny_regional_gwas_df.copy()]
        panel_labels = ["Only One"]  # Should have 2 labels

        with pytest.raises(ValueError, match="panel_labels"):
            regional_plotter.plot_stacked(
                gwas_dfs,
                chrom=1,
                start=1000000,
                end=2000000,
                panel_labels=panel_labels,
                show_recombination=False,
            )

    def test_stacked_mismatched_ld_reference_files_raises(
        self, regional_plotter, tiny_regional_gwas_df
    ):
        """Mismatched ld_reference_files length raises ValueError."""
        gwas_dfs = [tiny_regional_gwas_df, tiny_regional_gwas_df.copy()]
        ld_reference_files = ["/path/to/file1"]  # Only 1, but 2 gwas_dfs

        with pytest.raises(ValueError, match="ld_reference_files"):
            regional_plotter.plot_stacked(
                gwas_dfs,
                chrom=1,
                start=1000000,
                end=2000000,
                ld_reference_files=ld_reference_files,
                show_recombination=False,
            )


class TestFinemappingManyCredibleSets:
    """Test fine-mapping plot with many credible sets cycles colors correctly."""

    def test_finemapping_12_credible_sets_succeeds(
        self, regional_plotter, small_regional_gwas_df
    ):
        """Fine-mapping with >10 credible sets should cycle colors without error.

        The CREDIBLE_SET_COLORS palette has 10 colors, so 12 credible sets
        requires color cycling. This tests that modulo indexing works.
        """
        # Create fine-mapping data with 12 credible sets
        n_variants = 60  # 5 variants per credible set
        positions = list(range(1000000, 1000000 + n_variants * 10000, 10000))
        credible_sets = [((i // 5) % 12) + 1 for i in range(n_variants)]

        finemapping_df = pd.DataFrame(
            {
                "pos": positions,
                "pip": [0.8 if i % 5 == 0 else 0.1 for i in range(n_variants)],
                "cs": credible_sets,
            }
        )

        fig = regional_plotter.plot_stacked(
            [small_regional_gwas_df],
            chrom=1,
            start=900000,
            end=1700000,
            show_recombination=False,
            finemapping_df=finemapping_df,
            finemapping_cs_col="cs",
        )

        assert fig is not None

    def test_credible_set_color_cycling(self):
        """Verify get_credible_set_color cycles correctly for cs > 10."""
        from pylocuszoom.colors import CREDIBLE_SET_COLORS, get_credible_set_color

        # Test that colors cycle after 10
        for cs_id in range(1, 25):
            color = get_credible_set_color(cs_id)
            expected_idx = (cs_id - 1) % len(CREDIBLE_SET_COLORS)
            expected_color = CREDIBLE_SET_COLORS[expected_idx]
            assert color == expected_color, (
                f"CS {cs_id}: got {color}, expected {expected_color}"
            )


class TestRegionalPlotColumnValidation:
    """Test column validation in regional plots."""

    def test_plot_missing_pos_col_raises(self, regional_plotter):
        """Plot with missing position column should raise ValidationError."""
        from pylocuszoom.exceptions import ValidationError

        df = pd.DataFrame(
            {
                "rs": ["rs1", "rs2"],
                "p_wald": [1e-8, 0.01],
                # missing 'ps' column
            }
        )

        with pytest.raises(ValidationError, match="ps"):
            regional_plotter.plot(
                df,
                chrom=1,
                start=1000000,
                end=2000000,
            )

    def test_plot_missing_p_col_raises(self, regional_plotter):
        """Plot with missing p-value column should raise ValidationError."""
        from pylocuszoom.exceptions import ValidationError

        df = pd.DataFrame(
            {
                "rs": ["rs1", "rs2"],
                "ps": [1100000, 1900000],
                # missing 'p_wald' column
            }
        )

        with pytest.raises(ValidationError, match="p_wald"):
            regional_plotter.plot(
                df,
                chrom=1,
                start=1000000,
                end=2000000,
            )
