"""Tests for gene track visualization module."""

import matplotlib.pyplot as plt
import pandas as pd
import pytest
from hypothesis import given
from hypothesis import settings as hyp_settings

from pylocuszoom.colors import STRAND_COLORS
from pylocuszoom.config import RegionConfig
from pylocuszoom.gene_track import (
    assign_gene_positions,
    filter_genes_by_region,
    get_nearest_gene,
)
from pylocuszoom.panels.genes import GenePanel
from tests.strategies import gene_dataframes


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

    def test_complex_overlaps_no_same_row_collision(self):
        """Genes with complex overlaps should never share a row if they overlap.

        Regression test: the algorithm was incrementing row once per conflict
        without rechecking earlier rows, which could place overlapping genes
        in the same row.
        """
        # Create a scenario where naive single-pass would fail:
        # - Gene A: 100-300 (row 0)
        # - Gene B: 200-400 (row 1, conflicts with A)
        # - Gene C: 250-450 (row 2, conflicts with both A and B)
        # - Gene D: 150-350 (should be row 3, conflicts with all three)
        # BUT sorted by start: A, D, B, C - different order
        genes_df = pd.DataFrame(
            {
                "start": [100, 150, 200, 250],
                "end": [300, 350, 400, 450],
                "gene_name": ["GENE_A", "GENE_D", "GENE_B", "GENE_C"],
            }
        )
        # Must be sorted by start for algorithm
        genes_df = genes_df.sort_values("start")

        # Use a small region to ensure label buffer doesn't affect test
        positions = assign_gene_positions(genes_df, 0, 100000)

        # Verify no two overlapping genes share a row
        genes_with_rows = list(zip(genes_df["start"], genes_df["end"], positions))
        for i, (s1, e1, r1) in enumerate(genes_with_rows):
            for s2, e2, r2 in genes_with_rows[i + 1 :]:
                if r1 == r2:
                    # If same row, they shouldn't overlap
                    overlaps = not (e1 <= s2 or e2 <= s1)
                    assert not overlaps, (
                        f"Overlapping genes placed in same row {r1}: "
                        f"({s1}-{e1}) and ({s2}-{e2})"
                    )

    def test_three_overlapping_genes_correct_rows(self):
        """Three overlapping genes that would fail naive single-pass.

        Regression test: previous algorithm only incremented row counter
        without rechecking all occupied entries after increment.
        """
        # All three genes overlap each other significantly
        genes_df = pd.DataFrame(
            {
                "start": [1000, 1500, 2000],
                "end": [4000, 4500, 5000],
                "gene_name": ["GENE_A", "GENE_B", "GENE_C"],
            }
        )
        positions = assign_gene_positions(genes_df, 0, 100000)

        # Each gene should be in a different row
        assert len(set(positions)) == 3, (
            f"Three overlapping genes should be in 3 different rows, got {positions}"
        )
        # They should be in rows 0, 1, 2
        assert sorted(positions) == [0, 1, 2]


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

    def test_window_search_is_the_region_filter(self, genes_df):
        """The candidate set is the region filter centred on the position."""
        pos, window = 175000, 50000
        candidates = filter_genes_by_region(genes_df, 1, pos - window, pos + window)

        assert set(candidates["gene_name"]) == {"GENE_A", "GENE_B"}
        assert get_nearest_gene(genes_df, chrom=1, pos=pos, window=window) in set(
            candidates["gene_name"]
        )


def _draw(backend, ax, genes_df, chrom, start, end):
    """Draw a gene track for the region straight from a raw gene frame."""
    region = RegionConfig(chrom=str(chrom), start=start, end=end)
    GenePanel.from_genes(genes_df, region, None).draw(backend, ax)


class TestStrandColors:
    """Strand color constants are used by GenePanel.draw."""

    def test_strand_color_constants(self):
        assert "+" in STRAND_COLORS
        assert "-" in STRAND_COLORS
        assert None in STRAND_COLORS
        for color in STRAND_COLORS.values():
            assert color.startswith("#")
            assert len(color) == 7


# =============================================================================
# Property-Based Tests (Hypothesis)
# =============================================================================


class TestGeneTrackProperties:
    """Property-based tests for gene track rendering."""

    @given(gene_dataframes(min_genes=0, max_genes=30))
    def test_assign_positions_returns_valid_rows(self, genes_df):
        """Row assignments should always be non-negative integers."""
        if len(genes_df) == 0:
            positions = assign_gene_positions(genes_df, 0, 10000)
            assert positions == []
            return

        start = int(genes_df["start"].min())
        end = int(genes_df["end"].max())
        positions = assign_gene_positions(genes_df, start, end)

        assert len(positions) == len(genes_df)
        assert all(isinstance(p, int) and p >= 0 for p in positions)

    @given(gene_dataframes(min_genes=2, max_genes=20))
    def test_overlapping_genes_never_share_row(self, genes_df):
        """Overlapping genes should never be assigned to the same row."""
        if len(genes_df) < 2:
            return

        start = int(genes_df["start"].min())
        end = int(genes_df["end"].max())
        positions = assign_gene_positions(genes_df, start, end)

        # Check all pairs
        for i in range(len(genes_df)):
            for j in range(i + 1, len(genes_df)):
                gene_i_start = genes_df.iloc[i]["start"]
                gene_i_end = genes_df.iloc[i]["end"]
                gene_j_start = genes_df.iloc[j]["start"]
                gene_j_end = genes_df.iloc[j]["end"]

                # Check if genes overlap
                overlaps = gene_i_start < gene_j_end and gene_j_start < gene_i_end

                if overlaps:
                    assert positions[i] != positions[j], (
                        f"Overlapping genes at rows {i} and {j} share row {positions[i]}"
                    )

    @hyp_settings(max_examples=10, deadline=None)
    @given(gene_dataframes(min_genes=1, max_genes=15))
    def test_gene_panel_draws(self, genes_df):
        """GenePanel.draw should render without crashing."""
        if len(genes_df) == 0:
            return

        from pylocuszoom.backends.matplotlib_backend import MatplotlibBackend

        chrom = genes_df["chr"].iloc[0]
        start = int(genes_df["start"].min())
        end = int(genes_df["end"].max())

        fig, ax = plt.subplots()
        _draw(MatplotlibBackend(), ax, genes_df, chrom, start, end)


class TestStrandNaNHandling:
    """Regression: NaN strand should not silently render reversed arrows.

    Pre-fix bug: _draw_strand_arrows_matplotlib was called whenever the
    'strand' column existed; inside it `if strand == "+"` is False for
    NaN, falling through to the else branch and pointing arrows the
    wrong way. STRAND_COLORS.get(NaN) also missed the dict and used
    a fallback inconsistently.
    """

    def test_nan_strand_does_not_raise(self):
        """Plotting genes with NaN strand renders without error."""
        import numpy as np

        genes_df = pd.DataFrame(
            {
                "chr": [1, 1],
                "start": [1_000_000, 1_500_000],
                "end": [1_100_000, 1_600_000],
                "gene_name": ["KNOWN", "UNKNOWN_STRAND"],
                "strand": ["+", np.nan],
            }
        )
        fig, ax = plt.subplots()
        from pylocuszoom.backends.matplotlib_backend import MatplotlibBackend

        _draw(MatplotlibBackend(), ax, genes_df, 1, 1_000_000, 2_000_000)

    def test_invalid_strand_string_treated_as_missing(self):
        """A garbage strand value (not '+' or '-') is treated like NaN."""
        genes_df = pd.DataFrame(
            {
                "chr": [1],
                "start": [1_000_000],
                "end": [1_100_000],
                "gene_name": ["FOO"],
                "strand": ["?"],
            }
        )
        fig, ax = plt.subplots()
        from pylocuszoom.backends.matplotlib_backend import MatplotlibBackend

        _draw(MatplotlibBackend(), ax, genes_df, 1, 1_000_000, 2_000_000)
