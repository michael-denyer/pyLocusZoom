"""The regional LocusZoomPlotter: options, region selection and lead resolution."""

from unittest.mock import patch

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from pylocuszoom import ColumnConfig, DisplayConfig, LDConfig, PanelInputs
from pylocuszoom.backends.composition import LD_LEGEND_TITLE
from pylocuszoom.colors import LEAD_SNP_COLOR
from pylocuszoom.exceptions import ValidationError
from pylocuszoom.plotter import LocusZoomPlotter
from tests.figure_probes import PROBES

DISPLAY = DisplayConfig(show_recombination=False, snp_labels=False)


def _lead_marker_positions(fig, panel=0):
    """Return the x of every lead-SNP marker drawn on one panel."""
    return PROBES["matplotlib"].marker_x(fig, panel, color=LEAD_SNP_COLOR)


def points_and_limit(fig, backend):
    if backend == "matplotlib":
        return fig.axes[0].collections[0].get_offsets().tolist(), fig.axes[
            0
        ].get_ylim()[1]
    if backend == "plotly":
        return list(
            map(list, zip(fig.data[0].x, fig.data[0].y))
        ), fig.layout.yaxis.range[1]
    ax = fig.children[0]
    data = ax.renderers[0].data_source.data
    return list(map(list, zip(data["x"], data["y"]))), ax.y_range.end


@pytest.mark.parametrize("backend", ["matplotlib", "plotly", "bokeh"])
def test_region_selects_points_and_scale_before_drawing(backend):
    frame = pd.DataFrame(
        {
            "chr": [2, 1, 1],
            "pos": [150, 150, 9_000_000],
            "p_value": [1e-8, 0.01, 1e-200],
            "rs": ["wrong", "wanted", "far"],
        }
    )
    original = frame.copy(deep=True)
    plotter = LocusZoomPlotter(species=None, backend=backend, log_level=None)
    fig = plotter.plot(frame, chrom=1, start=100, end=200, display=DISPLAY)
    points, maximum = points_and_limit(fig, backend)
    assert points == [[150, 2.0]]
    assert maximum == pytest.approx(2.3)
    pd.testing.assert_frame_equal(frame, original)
    plt.close("all")


def test_ld_lead_comes_from_selected_chromosome():
    frame = pd.DataFrame(
        {
            "chr": [2, 1],
            "pos": [150, 150],
            "p_value": [1e-8, 0.01],
            "rs": ["wrong", "wanted"],
        }
    )
    leads = []

    def calculate(**kwargs):
        leads.append(kwargs["lead_snp"])
        return pd.DataFrame({"SNP": ["wanted"], "R2": [1.0]})

    with patch("pylocuszoom._ld_plotting.calculate_ld", side_effect=calculate):
        LocusZoomPlotter(species=None, log_level=None).plot(
            frame,
            chrom=1,
            start=100,
            end=200,
            display=DISPLAY,
            ld=LDConfig(lead_pos=150, ld_reference_file="reference"),
        )
    assert leads == ["wanted"]


@pytest.mark.parametrize("override, expected", [(None, 150), ([175], 175)])
def test_stacked_resolves_shared_lead_and_per_panel_override(override, expected):
    frame = pd.DataFrame({"pos": [150, 175], "p_value": [0.01, 1e-8]})
    fig = LocusZoomPlotter(species=None, backend="plotly", log_level=None).plot_stacked(
        [frame],
        chrom=1,
        start=100,
        end=200,
        display=DISPLAY,
        ld=LDConfig(lead_pos=150),
        lead_positions=override,
    )
    assert list(fig.data[-1].x) == [expected]


@pytest.mark.parametrize("lead", [0, -1])
def test_stacked_rejects_invalid_lead_positions(lead):
    frame = pd.DataFrame({"pos": [150], "p_value": [0.01]})
    with pytest.raises(ValueError, match="greater than or equal to 1"):
        LocusZoomPlotter(species=None, log_level=None).plot_stacked(
            [frame], chrom=1, start=100, end=200, display=DISPLAY, lead_positions=[lead]
        )


def test_stacked_resolves_columns_per_frame():
    canonical = pd.DataFrame({"pos": [150], "p_value": [0.01]})
    legacy = pd.DataFrame({"ps": [160], "p_wald": [0.001]})
    with pytest.warns(DeprecationWarning):
        fig = LocusZoomPlotter(
            species=None, backend="plotly", log_level=None
        ).plot_stacked(
            [canonical, legacy], chrom=1, start=100, end=200, display=DISPLAY
        )
    assert list(fig.data[0].x) == [150]
    assert list(fig.data[2].x) == [160]


def test_duplicate_index_does_not_make_lead_ambiguous():
    frame = pd.DataFrame({"pos": [150, 175], "p_value": [0.01, 1e-8]}, index=[0, 0])
    fig = LocusZoomPlotter(species=None, backend="plotly", log_level=None).plot(
        frame, chrom=1, start=100, end=200, display=DISPLAY, columns=ColumnConfig()
    )
    assert list(fig.data[-1].x) == [175]


@pytest.mark.parametrize("lead_positions", [None, [150]])
def test_duplicate_positions_use_one_strongest_lead_row(lead_positions):
    frame = pd.DataFrame(
        {
            "chr": [1, 1],
            "pos": [150, 150],
            "p_value": [0.1, 1e-8],
            "rs": ["weak", "strong"],
        }
    )
    observed_leads = []

    def calculate(**kwargs):
        observed_leads.append(kwargs["lead_snp"])
        return pd.DataFrame({"SNP": ["weak", "strong"], "R2": [0.2, 1.0]})

    with patch("pylocuszoom._ld_plotting.calculate_ld", side_effect=calculate):
        fig = LocusZoomPlotter(
            species=None, backend="plotly", log_level=None
        ).plot_stacked(
            [frame],
            chrom=1,
            start=100,
            end=200,
            display=DISPLAY,
            ld_reference_files=["reference"],
            lead_positions=lead_positions,
        )
    assert observed_leads == ["strong"]
    lead_trace = next(
        trace
        for trace in fig.data
        if trace.marker.symbol == "diamond" and trace.x[0] is not None
    )
    assert list(lead_trace.x) == [150]
    assert list(lead_trace.y) == [8.0]
    assert "strong" in str(lead_trace.customdata)


def test_requested_lead_missing_from_region_warns(warning_records):
    frame = pd.DataFrame({"pos": [150], "p_value": [0.1], "rs": ["other"]})
    fig = LocusZoomPlotter(species=None, backend="plotly", log_level=None).plot(
        frame,
        chrom=1,
        start=100,
        end=200,
        display=DISPLAY,
        ld=LDConfig(lead_pos=175, ld_reference_file="unused"),
    )
    assert list(fig.data[0].x) == [150]
    assert any("175 not found in region" in message for message in warning_records)


def test_same_position_nonlead_variant_is_excluded_from_lead_labels():
    frame = pd.DataFrame(
        {"pos": [150, 150], "p_value": [0.1, 1e-8], "rs": ["weak", "strong"]}
    )
    fig = LocusZoomPlotter(species=None, log_level=None).plot(
        frame,
        chrom=1,
        start=100,
        end=200,
        display=DisplayConfig(show_recombination=False, label_top_n=2),
    )
    assert [text.get_text() for text in fig.axes[0].texts] == ["strong"]


class TestLocusZoomPlotterInit:
    """Tests for LocusZoomPlotter initialization."""

    def test_default_species_is_canine(self):
        """Default species should be canine."""
        plotter = LocusZoomPlotter()
        assert plotter.species.key == "canine"

    def test_custom_species(self):
        """Should accept custom species."""
        plotter = LocusZoomPlotter(species="feline")
        assert plotter.species.key == "feline"

    def test_custom_plink_path(self):
        """Should accept custom PLINK path."""
        plotter = LocusZoomPlotter(plink_path="/custom/plink")
        assert plotter.plink_path == "/custom/plink"

    def test_custom_threshold(self):
        """Should accept custom genomewide threshold."""
        plotter = LocusZoomPlotter(genomewide_threshold=5e-8)
        assert plotter.genomewide_threshold == 5e-8


class TestLocusZoomPlotterPlot:
    """Tests for LocusZoomPlotter.plot() method."""

    def test_recombination_overlay_is_on_by_default(
        self, canine_plotter, regional_gwas_df
    ):
        """A species with maps gets a recombination-rate axis over the association."""
        fig = canine_plotter.plot(
            regional_gwas_df,
            chrom=1,
            start=1000000,
            end=2000000,
        )
        assert [ax.get_ylabel() for ax in fig.axes] == [
            r"$-\log_{10}$ P",
            "Recombination rate (cM/Mb)",
        ]
        # The fixture's own map: its one in-region row, at its canfam3.1 position.
        (recomb_line,) = fig.axes[1].get_lines()
        assert list(recomb_line.get_xdata()) == [1500000]
        assert list(recomb_line.get_ydata()) == [1.0]

    def test_plots_with_gene_track(
        self, canine_plotter, regional_gwas_df, sample_genes_df
    ):
        """Should create plot with gene track when genes_df provided."""
        fig = canine_plotter.plot(
            regional_gwas_df,
            chrom=1,
            start=1000000,
            end=2000000,
            panels=PanelInputs(genes_df=sample_genes_df),
        )
        # Should have 2 axes (association + gene track)
        assert len(fig.axes) >= 2

    def test_highlights_lead_snp(self, canine_plotter, regional_gwas_df):
        """Should highlight lead SNP when lead_pos provided."""
        lead_pos = regional_gwas_df["pos"].iloc[0]
        fig = canine_plotter.plot(
            regional_gwas_df,
            chrom=1,
            start=1000000,
            end=2000000,
            ld=LDConfig(lead_pos=lead_pos),
        )
        assert _lead_marker_positions(fig) == [lead_pos]

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
            columns=ColumnConfig(pos_col="position", p_col="pvalue", rs_col="snp_id"),
            display=DisplayConfig(show_recombination=False),
        )
        assert set(PROBES["matplotlib"].marker_x(fig)) == {
            1100000.0,
            1500000.0,
            1900000.0,
        }
        assert _lead_marker_positions(fig) == [1100000.0]

    def test_with_precomputed_ld(self, canine_plotter, regional_gwas_df):
        """Should use pre-computed LD column when provided."""
        df = regional_gwas_df.copy()
        df["R2"] = np.random.default_rng(0).uniform(0, 1, len(df))

        fig = canine_plotter.plot(
            df, chrom=1, start=1000000, end=2000000, ld=LDConfig(ld_col="R2")
        )
        legend = fig.axes[0].get_legend()
        assert legend is not None
        assert legend.get_title().get_text() == LD_LEGEND_TITLE

    def test_with_recombination_data(
        self, canine_plotter, regional_gwas_df, sample_recomb_df
    ):
        """Should plot with recombination overlay when provided."""
        fig = canine_plotter.plot(
            regional_gwas_df,
            chrom=1,
            start=1000000,
            end=2000000,
            panels=PanelInputs(recomb_df=sample_recomb_df),
        )
        (recomb_line,) = fig.axes[1].get_lines()
        assert list(recomb_line.get_xdata()) == list(sample_recomb_df["pos"])
        assert list(recomb_line.get_ydata()) == list(sample_recomb_df["rate"])

    def test_disables_snp_labels(self, canine_plotter, regional_gwas_df):
        """Should not add labels when snp_labels=False."""
        fig = canine_plotter.plot(
            regional_gwas_df,
            chrom=1,
            start=1000000,
            end=2000000,
            display=DisplayConfig(snp_labels=False, show_recombination=False),
        )
        assert list(fig.axes[0].texts) == []

    def test_disables_recombination(self, canine_plotter, regional_gwas_df):
        """Should not show recombination when show_recombination=False."""
        fig = canine_plotter.plot(
            regional_gwas_df,
            chrom=1,
            start=1000000,
            end=2000000,
            display=DisplayConfig(show_recombination=False),
        )
        assert len(fig.axes) == 1


class TestFloatChromosomeColumn:
    """A float chromosome column filters like its integer twin."""

    @pytest.mark.parametrize("chroms", [[1, 1, 1], [1.0, 1.0, 1.0]])
    def test_regional_plot_draws_the_region_points(self, chroms):
        df = pd.DataFrame(
            {
                "chr": chroms,
                "pos": [1100, 1500, 1900],
                "p_value": [1e-3, 1e-9, 1e-5],
                "rs": ["a", "b", "c"],
            }
        )
        plotter = LocusZoomPlotter(species="canine", log_level=None)

        fig = plotter.plot(
            df,
            chrom=1,
            start=1000,
            end=2000,
            display=DisplayConfig(show_recombination=False),
        )

        assert set(PROBES["matplotlib"].marker_x(fig)) == {1100.0, 1500.0, 1900.0}


class TestPlotStackedEdgeCases:
    """Tests for plot_stacked() edge cases and error handling."""

    def test_plot_stacked_validates_eqtl_columns(
        self, canine_plotter, tiny_regional_gwas_df
    ):
        """plot_stacked() validates the eQTL frame instead of raising KeyError."""
        from pylocuszoom.eqtl import EQTLValidationError

        # eQTL data with wrong column names
        bad_eqtl_df = pd.DataFrame(
            {
                "position": [1500000],  # Should be 'pos'
                "pval": [1e-6],  # Should be 'p_value'
            }
        )

        with pytest.raises(EQTLValidationError):
            canine_plotter.plot_stacked(
                [tiny_regional_gwas_df],
                chrom=1,
                start=1000000,
                end=2000000,
                display=DisplayConfig(show_recombination=False),
                panels=PanelInputs(eqtl_df=bad_eqtl_df),
            )


class TestStackedPlotMismatchedLengths:
    """A per-panel list whose length differs from gwas_dfs raises ValueError."""

    @pytest.mark.parametrize(
        ("option", "value"),
        [
            ("lead_positions", [1500000, 1500000]),
            ("panel_labels", ["A", "B"]),
            ("ld_reference_files", ["/path/to/file1", "/path/to/file2"]),
        ],
    )
    def test_stacked_list_of_the_wrong_length_raises(
        self, regional_plotter, tiny_regional_gwas_df, option, value
    ):
        gwas_dfs = [tiny_regional_gwas_df] * 3

        with pytest.raises(ValueError, match=option):
            regional_plotter.plot_stacked(
                gwas_dfs,
                chrom=1,
                start=1000000,
                end=2000000,
                display=DisplayConfig(show_recombination=False),
                **{option: value},
            )


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
                "pos": [1_200_000, 1_800_000, 1_500_000, 1_900_000],
                "p_value": [1e-5, 1e-3, 1e-12, 1e-10],
            }
        )

        fig = canine_plotter.plot_stacked(
            [gwas_df],
            chrom=1,
            start=1_000_000,
            end=2_000_000,
            columns=ColumnConfig(pos_col="pos", p_col="p_value"),
            display=DisplayConfig(show_recombination=False),
        )

        assert _lead_marker_positions(fig) == [1_200_000], (
            "Lead must be chr1's strongest hit (1_200_000), not chr2's (1_500_000)"
        )


class TestLeadAutoDetectionAgreesAcrossEntryPoints:
    """With no lead given, plot() and plot_stacked([df]) pick the same lead:
    the strongest in-region p-value."""

    def test_plot_auto_detects_the_strongest_hit(self, canine_plotter):
        gwas_df = pd.DataFrame(
            {
                "rs": ["rs1", "rs2", "rs3"],
                "chr": [1, 1, 1],
                "pos": [1_200_000, 1_500_000, 1_800_000],
                "p_value": [1e-3, 1e-9, 1e-5],
            }
        )
        kwargs = dict(
            chrom=1,
            start=1_000_000,
            end=2_000_000,
            display=DisplayConfig(show_recombination=False),
        )
        single = canine_plotter.plot(gwas_df, **kwargs)
        stacked = canine_plotter.plot_stacked([gwas_df], **kwargs)

        assert _lead_marker_positions(single) == [1_500_000]
        assert _lead_marker_positions(stacked) == [1_500_000]


class TestLeadPosBoundary:
    """Pins the lead_pos=1 boundary.

    Smallest valid position (``ge=1`` in config) reaches the association
    panel intact, and the public API rejects ``0`` (1-based genomic coords).
    """

    def test_lead_pos_one_reaches_plot_association(self):
        """lead_pos=1 (smallest valid position) marks the SNP at position 1.

        The SNP at 1 is not the strongest, so a falsy check that dropped the
        lead to None would move the marker to the auto-detected hit.
        """
        plotter = LocusZoomPlotter(species="canine", log_level=None)
        gwas_df = pd.DataFrame(
            {
                "rs": ["rs_lead", "rs2", "rs3"],
                "pos": [1, 100_000, 200_000],
                "p_value": [1e-3, 1e-8, 1e-5],
            }
        )

        fig = plotter.plot(
            gwas_df,
            chrom=1,
            start=1,
            end=300_000,
            columns=ColumnConfig(pos_col="pos", p_col="p_value"),
            display=DisplayConfig(show_recombination=False),
            ld=LDConfig(lead_pos=1),
        )

        assert _lead_marker_positions(fig) == [1.0]

    def test_lead_pos_zero_rejected_at_api(self):
        """Public API enforces genomic coords are 1-based; lead_pos=0 rejected."""
        plotter = LocusZoomPlotter(species="canine", log_level=None)
        gwas_df = pd.DataFrame(
            {
                "rs": ["rs1", "rs2"],
                "pos": [100_000, 200_000],
                "p_value": [1e-8, 1e-5],
            }
        )

        with pytest.raises(ValidationError, match="greater than or equal to 1"):
            plotter.plot(
                gwas_df,
                chrom=1,
                start=1,
                end=300_000,
                columns=ColumnConfig(pos_col="pos", p_col="p_value"),
                display=DisplayConfig(show_recombination=False),
                ld=LDConfig(lead_pos=0),
            )
