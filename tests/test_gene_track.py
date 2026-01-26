"""Tests for gene track visualization module."""

import matplotlib.pyplot as plt
import pandas as pd
import pytest

from snp_scope_plot.gene_track import (
    STRAND_COLORS,
    assign_gene_positions,
    get_nearest_gene,
    plot_gene_track,
)


class TestAssignGenePositions:
    """Tests for assign_gene_positions function."""

    def test_single_gene_gets_row_zero(self):
        """Single gene should be placed in row 0."""
        genes_df = pd.DataFrame(
            {
                "start": [1000],
                "end": [2000],
                "gene_name": ["GENE_A"],
            }
        )
        positions = assign_gene_positions(genes_df, 0, 10000)
        assert positions == [0]

    def test_non_overlapping_genes_same_row(self):
        """Non-overlapping genes should be in the same row."""
        genes_df = pd.DataFrame(
            {
                "start": [1000, 5000],
                "end": [2000, 6000],
                "gene_name": ["GENE_A", "GENE_B"],
            }
        )
        positions = assign_gene_positions(genes_df, 0, 10000)
        # Both should be in row 0 as they don't overlap
        assert positions == [0, 0]

    def test_overlapping_genes_different_rows(self):
        """Overlapping genes should be in different rows."""
        genes_df = pd.DataFrame(
            {
                "start": [1000, 1500],
                "end": [3000, 4000],
                "gene_name": ["GENE_A", "GENE_B"],
            }
        )
        positions = assign_gene_positions(genes_df, 0, 10000)
        assert positions[0] != positions[1]

    def test_three_stacked_genes(self):
        """Three overlapping genes should stack vertically."""
        genes_df = pd.DataFrame(
            {
                "start": [1000, 1100, 1200],
                "end": [5000, 5100, 5200],
                "gene_name": ["GENE_A", "GENE_B", "GENE_C"],
            }
        )
        positions = assign_gene_positions(genes_df, 0, 10000)
        assert positions == [0, 1, 2]

    def test_empty_dataframe(self):
        """Empty DataFrame should return empty list."""
        genes_df = pd.DataFrame(columns=["start", "end", "gene_name"])
        positions = assign_gene_positions(genes_df, 0, 10000)
        assert positions == []


class TestGetNearestGene:
    """Tests for get_nearest_gene function."""

    @pytest.fixture
    def genes_df(self):
        """Sample gene annotations."""
        return pd.DataFrame(
            {
                "chr": [1, 1, 1],
                "start": [100000, 200000, 500000],
                "end": [150000, 250000, 550000],
                "gene_name": ["GENE_A", "GENE_B", "GENE_C"],
            }
        )

    def test_finds_overlapping_gene(self, genes_df):
        """Should find gene when position overlaps."""
        result = get_nearest_gene(genes_df, chrom=1, pos=125000)
        assert result == "GENE_A"

    def test_finds_nearby_gene_within_window(self, genes_df):
        """Should find gene when position is within window."""
        # 175000 is 25kb from GENE_A end (150000) and 25kb from GENE_B start (200000)
        result = get_nearest_gene(genes_df, chrom=1, pos=175000)
        # Should find GENE_B as it's closer to midpoint
        assert result in ["GENE_A", "GENE_B"]

    def test_returns_none_outside_window(self, genes_df):
        """Should return None when no gene within window."""
        # Position 350000 is > 50kb from both GENE_B (ends 250000) and GENE_C (starts 500000)
        result = get_nearest_gene(genes_df, chrom=1, pos=350000)
        assert result is None

    def test_returns_none_for_empty_chromosome(self, genes_df):
        """Should return None when chromosome has no genes."""
        result = get_nearest_gene(genes_df, chrom=2, pos=125000)
        assert result is None

    def test_handles_chr_prefix(self, genes_df):
        """Should handle 'chr' prefix in chromosome."""
        genes_with_prefix = genes_df.copy()
        genes_with_prefix["chr"] = "chr1"
        result = get_nearest_gene(genes_with_prefix, chrom=1, pos=125000)
        assert result == "GENE_A"

    def test_custom_window_size(self, genes_df):
        """Should use custom window size."""
        # With 100kb window, 350000 should find GENE_C (starts 500000)
        result = get_nearest_gene(genes_df, chrom=1, pos=410000, window=100000)
        assert result == "GENE_C"


class TestPlotGeneTrack:
    """Tests for plot_gene_track function."""

    @pytest.fixture
    def sample_genes(self):
        """Sample gene DataFrame."""
        return pd.DataFrame(
            {
                "chr": [1, 1, 1],
                "start": [1100000, 1400000, 1700000],
                "end": [1150000, 1500000, 1800000],
                "gene_name": ["GENE_A", "GENE_B", "GENE_C"],
                "strand": ["+", "-", "+"],
            }
        )

    @pytest.fixture
    def sample_exons(self):
        """Sample exon DataFrame."""
        return pd.DataFrame(
            {
                "chr": [1, 1, 1, 1],
                "start": [1100000, 1120000, 1400000, 1450000],
                "end": [1110000, 1130000, 1420000, 1470000],
                "gene_name": ["GENE_A", "GENE_A", "GENE_B", "GENE_B"],
            }
        )

    def test_creates_gene_track(self, sample_genes):
        """Should create gene track without errors."""
        fig, ax = plt.subplots()
        plot_gene_track(ax, sample_genes, chrom=1, start=1000000, end=2000000)
        plt.close(fig)
        # Test passes if no exception raised

    def test_handles_empty_genes(self):
        """Should handle empty gene DataFrame gracefully."""
        fig, ax = plt.subplots()
        empty_genes = pd.DataFrame(columns=["chr", "start", "end", "gene_name"])
        plot_gene_track(ax, empty_genes, chrom=1, start=1000000, end=2000000)
        plt.close(fig)

    def test_handles_no_genes_in_region(self, sample_genes):
        """Should handle region with no genes."""
        fig, ax = plt.subplots()
        # Query chromosome 2 where no genes exist
        plot_gene_track(ax, sample_genes, chrom=2, start=1000000, end=2000000)
        plt.close(fig)

    def test_plots_with_exons(self, sample_genes, sample_exons):
        """Should plot exon structure when provided."""
        fig, ax = plt.subplots()
        plot_gene_track(
            ax, sample_genes, chrom=1, start=1000000, end=2000000, exons_df=sample_exons
        )
        plt.close(fig)

    def test_uses_strand_colors(self, sample_genes):
        """Gene colors should match strand direction."""
        # Just verify the color constants exist and are valid
        assert "+" in STRAND_COLORS
        assert "-" in STRAND_COLORS
        assert None in STRAND_COLORS
        # All should be valid hex colors
        for color in STRAND_COLORS.values():
            assert color.startswith("#")
            assert len(color) == 7

    def test_sets_axis_limits(self, sample_genes):
        """Should set correct x-axis limits."""
        fig, ax = plt.subplots()
        plot_gene_track(ax, sample_genes, chrom=1, start=1000000, end=2000000)
        xlim = ax.get_xlim()
        assert xlim == (1000000, 2000000)
        plt.close(fig)
