"""Automatic gene fetching for the regional plot's gene track."""

from unittest.mock import Mock, patch

import pandas as pd
import pytest

from pylocuszoom import DisplayConfig, PanelInputs, ValidationError
from pylocuszoom._gene_source import GeneAnnotations
from pylocuszoom.plotter import LocusZoomPlotter
from tests.reference_mocks import (
    gene_transcript_exon_payload,
    ok_response,
    refseq_payload,
)


class TestSuppliedExons:
    def test_malformed_exons_raise_naming_the_frame(
        self, small_regional_gwas_df, sample_genes_df
    ):
        """An exon frame with 'chrom' for 'chr' used to raise a bare KeyError."""
        exons = pd.DataFrame(
            {
                "chrom": [1],
                "start": [1110000],
                "end": [1120000],
                "gene_name": ["GENE_A"],
            }
        )

        with pytest.raises(ValidationError, match="(?s)exons_df.*'chr'"):
            LocusZoomPlotter(species=None, log_level=None).plot(
                small_regional_gwas_df,
                chrom=1,
                start=1000000,
                end=2000000,
                display=DisplayConfig(show_recombination=False),
                panels=PanelInputs(genes_df=sample_genes_df, exons_df=exons),
            )


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
        """Label the gene track with the genes the reference source returned."""
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

        gene_track_ax = fig.get_axes()[1]
        assert [text.get_text() for text in gene_track_ax.texts] == ["GENE1", "GENE2"]

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
                display=DisplayConfig(show_recombination=False),
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
            pytest.warns(
                UserWarning, match=r"chr1:1000000-2000000.*UCSC.*503"
            ) as caught,
        ):
            fig = plotter.plot(
                small_regional_gwas_df,
                chrom=1,
                start=1000000,
                end=2000000,
                display=DisplayConfig(show_recombination=False),
            )

        assert [w.filename for w in caught] == [__file__]
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
                display=DisplayConfig(show_recombination=False, auto_genes=True),
            )

        mock_fetch.assert_called_once()
        assert len(fig.get_axes()) == 2

    def test_plot_auto_genes_disabled_by_default(self, small_regional_gwas_df):
        """Without auto_genes a plot reaches no reference source."""
        plotter = LocusZoomPlotter(species="canine", log_level=None)

        with patch("pylocuszoom.plotter.get_genes_for_build") as mock_fetch:
            plotter.plot(
                small_regional_gwas_df,
                chrom=1,
                start=1000000,
                end=2000000,
                display=DisplayConfig(show_recombination=False),
            )

        assert not mock_fetch.called, "auto_genes is off, so no gene fetch is allowed"

    def test_plot_auto_genes_respects_explicit_genes_df(
        self, small_regional_gwas_df, two_gene_track_df
    ):
        """Draw the caller's own genes, and fetch nothing, when genes_df is given."""
        plotter = LocusZoomPlotter(species="human", log_level=None, auto_genes=True)

        with patch("pylocuszoom.plotter.get_genes_for_build") as mock_fetch:
            fig = plotter.plot(
                small_regional_gwas_df,
                chrom=1,
                start=1000000,
                end=2000000,
                panels=PanelInputs(genes_df=two_gene_track_df),
            )

            mock_fetch.assert_not_called()

        gene_track_ax = fig.get_axes()[1]
        assert [
            (text.get_text(), text.get_position()[0]) for text in gene_track_ax.texts
        ] == [("GENE1", 1150000.0), ("GENE2", 1350000.0)]
