"""Liftover correctness against pyliftover's real coordinate contract.

pyliftover takes and returns 0-based positions, while every position column in
this package is 1-based. The chains below put a block boundary exactly where an
off-by-one query lands in the wrong block, so these tests use real pyliftover
rather than a fake that could encode the wrong contract.
"""

import gzip
import warnings

import pandas as pd
import pytest
from pyliftover import LiftOver

from pylocuszoom import (
    ColumnConfig,
    CoordinateLifter,
    DisplayConfig,
    LDConfig,
    LiftoverConfig,
    LocusZoomPlotter,
    liftover_region,
)
from pylocuszoom._liftover import InMemoryLifter, chain_lifter
from pylocuszoom.exceptions import (
    DataDownloadError,
    OptionalDependencyMissing,
    ValidationError,
)
from pylocuszoom.genome_build import GENOME_BUILDS, GenomeBuild
from pylocuszoom.recombination import get_recombination_rate_for_region

# chr1 [0, 1000) -> chr1 [100, 1100); chr1 [1000, 2000) -> chr5 [0, 1000).
# 1-based chr1:1000 is 0-based 999, the last base of the first block.
SPLIT_CHAIN = (
    "chain 1000 chr1 5000 + 0 1000 chr1 5000 + 100 1100 1\n1000\n\n"
    "chain 1000 chr1 5000 + 1000 2000 chr5 5000 + 0 1000 2\n1000\n\n"
)


@pytest.fixture
def split_lifter(tmp_path):
    path = tmp_path / "split.chain"
    path.write_text(SPLIT_CHAIN)
    return LiftOver(str(path))


class TestLiftoverRegion:
    """The public region liftover behind regional plots across builds."""

    @pytest.fixture
    def region_df(self):
        return pd.DataFrame(
            {
                "pos": [100, 200, 300],
                "p_value": [1e-9, 1e-3, 0.2],
                "rs": ["rsA", "rsB", "rsC"],
                "R2": [1.0, 0.5, 0.1],
            }
        )

    def test_real_chain_lifts_block_edge_and_drops_other_chromosome(self, split_lifter):
        df = pd.DataFrame({"pos": [1000, 1500], "p_value": [1e-8, 1e-3]})

        result = liftover_region(df, chrom=1, lifter=split_lifter, lead_pos=1000)

        assert result.lifted_df["pos"].tolist() == [1100]
        assert result.lead_pos == 1100
        assert (result.n_lifted, result.n_cross_chrom, result.n_dropped) == (1, 1, 1)

    def test_pyliftover_satisfies_the_protocol(self, split_lifter):
        assert isinstance(split_lifter, CoordinateLifter)

    def test_keeps_other_columns_and_row_order(self, region_df):
        lifter = InMemoryLifter(
            {("chr1", 99): 1099, ("chr1", 199): 1199, ("chr1", 299): 1299}
        )

        result = liftover_region(region_df, chrom=1, lifter=lifter, lead_pos=100)

        assert result.lifted_df["pos"].tolist() == [1100, 1200, 1300]
        assert result.lifted_df["R2"].tolist() == [1.0, 0.5, 0.1]
        assert (result.start, result.end, result.lead_pos) == (1100, 1300, 1100)
        assert result.is_collinear is True
        assert result.source_lead_pos == 100

    def test_counts_every_drop_reason(self, region_df):
        df = pd.concat([region_df, region_df.iloc[[0]].assign(pos=400)])
        lifter = InMemoryLifter(
            {
                ("chr1", 99): 1099,
                ("chr1", 199): [1199, 4999],
                ("chr1", 299): ("chr5", 1299),
            }
        )

        result = liftover_region(df, chrom="chr1", lifter=lifter)

        assert result.lifted_df["rs"].tolist() == ["rsA"]
        assert (result.n_unmapped, result.n_multimapped, result.n_cross_chrom) == (
            1,
            1,
            1,
        )
        assert (result.n_dropped, result.n_input) == (3, 4)

    def test_detects_local_rearrangement(self, region_df):
        lifter = InMemoryLifter(
            {("chr1", 99): 2999, ("chr1", 199): 1999, ("chr1", 299): 999}
        )

        result = liftover_region(region_df, chrom=1, lifter=lifter)

        assert result.is_collinear is False

    def test_nothing_lifts(self, region_df):
        result = liftover_region(
            region_df, chrom=1, lifter=InMemoryLifter({("chr1", 5): 5})
        )

        assert result.lifted_df.empty
        assert (result.start, result.end) == (None, None)
        assert result.n_unmapped == 3

    def test_reports_a_chromosome_the_chain_lacks_without_warning(
        self, region_df, recwarn
    ):
        result = liftover_region(
            region_df, chrom=1, lifter=InMemoryLifter({("chr2", 0): 0})
        )

        assert result.chain_has_chrom is False
        assert result.lifted_df.empty
        assert list(recwarn) == []

    @pytest.mark.parametrize("chrom", [39, "39", 41, "XY", "X", "chrX"])
    def test_canine_plink_x_codes_query_chr_x(self, chrom, tmp_path):
        path = tmp_path / "x.chain"
        path.write_text(
            "chain 1000 chrX 5000 + 0 1000 chrX 5000 + 100 1100 1\n1000\n\n"
        )
        df = pd.DataFrame({"pos": [1000], "p_value": [1e-9]})

        result = liftover_region(
            df,
            chrom=chrom,
            lifter=LiftOver(str(path)),
            build="canfam3",
        )

        assert result.lifted_df["pos"].tolist() == [1100]

    @pytest.mark.parametrize("chrom", [19, 21, "XY"])
    def test_feline_plink_x_codes_query_chr_x(self, chrom):
        lifter = InMemoryLifter({("chrX", 99): 1099})
        df = pd.DataFrame({"pos": [100], "p_value": [1e-9]})

        result = liftover_region(df, chrom=chrom, lifter=lifter, build="felCat9")

        assert result.lifted_df["pos"].tolist() == [1100]

    def test_numeric_x_code_is_not_aliased_without_a_build(self):
        lifter = InMemoryLifter({("chrX", 99): 1099, ("chr39", 0): 0})
        df = pd.DataFrame({"pos": [100], "p_value": [1e-9]})

        result = liftover_region(df, chrom=39, lifter=lifter)

        assert result.lifted_df.empty


class TestChainLifter:
    """Registered chains download once into the cache and load through load_chain."""

    CANFAM3, CANFAM4 = GENOME_BUILDS["canfam3"], GENOME_BUILDS["canfam4"]

    @pytest.fixture
    def downloads(self, monkeypatch):
        """Record every chain download, writing SPLIT_CHAIN to the destination."""
        urls = []

        def download(url, dest, desc=None):
            urls.append(url)
            dest.write_bytes(gzip.compress(SPLIT_CHAIN.encode()))

        monkeypatch.setattr("pylocuszoom._liftover.download_file", download)
        return urls

    def test_downloads_a_missing_chain_into_the_cache(self, cache_home, downloads):
        lifter = chain_lifter(self.CANFAM3, self.CANFAM4)

        assert downloads == [self.CANFAM3.chain_url(self.CANFAM4)]
        assert (cache_home / "liftover" / "canFam3ToCanFam4.over.chain.gz").exists()
        assert lifter.convert_coordinate("chr1", 999)[0][:2] == ("chr1", 1099)

    def test_reuses_a_cached_chain(self, cache_home, downloads):
        (cache_home / "liftover").mkdir(parents=True)
        (cache_home / "liftover" / "canFam3ToCanFam4.over.chain.gz").write_bytes(
            gzip.compress(SPLIT_CHAIN.encode())
        )

        chain_lifter(self.CANFAM3, self.CANFAM4)

        assert downloads == []

    def test_each_registered_chain_downloads_its_own_url(self, cache_home, downloads):
        url = "https://example.org/chains/fooToBar.over.chain"
        source = GenomeBuild(
            key="foo", species="x", assembly_name="Foo", liftover_chains=(("bar", url),)
        )
        target = GenomeBuild(key="bar", species="x", assembly_name="Bar")

        chain_lifter(source, target)

        assert downloads == [url]

    def test_an_unregistered_pair_is_a_validation_error(self, cache_home, downloads):
        with pytest.raises(ValidationError, match="No liftover chain"):
            chain_lifter(self.CANFAM4, self.CANFAM3)
        assert downloads == []

    def test_a_chain_that_stays_unreadable_is_removed_and_reported(
        self, cache_home, monkeypatch
    ):
        def download(url, dest, desc=None):
            dest.write_bytes(b"<html>502</html>")

        monkeypatch.setattr("pylocuszoom._liftover.download_file", download)

        with pytest.raises(DataDownloadError, match="unreadable"):
            chain_lifter(self.CANFAM3, self.CANFAM4)
        assert not (cache_home / "liftover" / "canFam3ToCanFam4.over.chain.gz").exists()

    def test_missing_pyliftover_is_a_typed_error(self, tmp_path, monkeypatch):
        from pylocuszoom._liftover import load_chain

        monkeypatch.setitem(__import__("sys").modules, "pyliftover", None)
        path = tmp_path / "split.chain"
        path.write_text(SPLIT_CHAIN)

        with pytest.raises(OptionalDependencyMissing, match="pip install pyliftover"):
            load_chain(path)


class TestLiftoverConfig:
    def test_rejects_lifter_and_chain_together(self):
        with pytest.raises(ValueError, match="either lifter or chain_path"):
            LiftoverConfig(lifter=InMemoryLifter({}), chain_path="x.chain")

    def test_lift_recombination_needs_a_chain(self):
        with pytest.raises(ValueError, match="lift_recombination requires"):
            LiftoverConfig(lift_recombination=True)

    def test_rejects_an_object_that_is_not_a_lifter(self):
        with pytest.raises(ValueError):
            LiftoverConfig(lifter=object())

    def test_chain_path_loads_pyliftover(self, tmp_path):
        path = tmp_path / "split.chain"
        path.write_text(SPLIT_CHAIN)

        lifter = LiftoverConfig(chain_path=path).resolve()

        assert lifter.convert_coordinate("chr1", 999)[0][:2] == ("chr1", 1099)

    def test_default_is_no_liftover(self):
        assert LiftoverConfig().resolve() is None


class TestPlotAcrossBuilds:
    """LocusZoomPlotter.plot lifts source-build sumstats before drawing."""

    @pytest.fixture
    def plotter(self):
        return LocusZoomPlotter(species="canine", genome_build="canfam4")

    @pytest.fixture
    def source_gwas_df(self):
        return pd.DataFrame(
            {
                "chr": [1, 1, 1],
                "ps": [1_000, 2_000, 3_000],
                "p_wald": [1e-9, 1e-3, 0.2],
                "rs": ["rsA", "rsB", "rsC"],
            }
        )

    LIFTER = InMemoryLifter(
        {("chr1", 999): 10_999, ("chr1", 1_999): 11_999, ("chr1", 2_999): 12_999}
    )

    def test_plots_lifted_positions_with_requested_margins(
        self, plotter, source_gwas_df
    ):
        fig = plotter.plot(
            source_gwas_df,
            chrom=1,
            start=500,
            end=3_500,
            columns=ColumnConfig(pos_col="ps", p_col="p_wald"),
            ld=LDConfig(lead_pos=1_000),
            display=DisplayConfig(show_recombination=False, snp_labels=False),
            liftover=LiftoverConfig(lifter=self.LIFTER),
        )

        ax = fig.axes[0]
        xs = {x for coll in ax.collections for x, _ in coll.get_offsets()}
        assert xs == {11_000, 12_000, 13_000}
        assert ax.get_xlim() == (10_500, 13_500)

    @pytest.mark.parametrize("chrom_col", ["chrom", None])
    def test_lifts_a_frame_by_its_configured_chromosome_column(
        self, plotter, source_gwas_df, chrom_col
    ):
        frame = source_gwas_df.drop(columns="chr")
        if chrom_col is not None:
            frame[chrom_col] = source_gwas_df["chr"]
        fig = plotter.plot(
            frame,
            chrom=1,
            start=500,
            end=3_500,
            columns=ColumnConfig(chrom_col=chrom_col, pos_col="ps", p_col="p_wald"),
            display=DisplayConfig(show_recombination=False, snp_labels=False),
            liftover=LiftoverConfig(lifter=self.LIFTER),
        )

        xs = {x for coll in fig.axes[0].collections for x, _ in coll.get_offsets()}
        assert xs == {11_000, 12_000, 13_000}

    def test_raises_when_nothing_lifts(self, plotter, source_gwas_df):
        with pytest.raises(ValueError, match="No SNP in chr1:500-3500 lifted"):
            plotter.plot(
                source_gwas_df,
                chrom=1,
                start=500,
                end=3_500,
                columns=ColumnConfig(pos_col="ps", p_col="p_wald"),
                display=DisplayConfig(show_recombination=False),
                liftover=LiftoverConfig(lifter=InMemoryLifter({("chr1", 0): 0})),
            )

    def test_warns_when_lead_does_not_lift(self, plotter, source_gwas_df):
        lifter = InMemoryLifter({("chr1", 1_999): 11_999, ("chr1", 2_999): 12_999})
        with pytest.warns(UserWarning, match="Lead SNP at chr1:1000 did not lift"):
            plotter.plot(
                source_gwas_df,
                chrom=1,
                start=500,
                end=3_500,
                columns=ColumnConfig(pos_col="ps", p_col="p_wald"),
                ld=LDConfig(lead_pos=1_000),
                display=DisplayConfig(show_recombination=False, snp_labels=False),
                liftover=LiftoverConfig(lifter=lifter),
            )

    def test_warns_when_the_region_is_rearranged(self, plotter, source_gwas_df):
        lifter = InMemoryLifter(
            {("chr1", 999): 12_999, ("chr1", 1_999): 11_999, ("chr1", 2_999): 10_999}
        )
        with pytest.warns(UserWarning, match="rearranged between builds"):
            plotter.plot(
                source_gwas_df,
                chrom=1,
                start=500,
                end=3_500,
                columns=ColumnConfig(pos_col="ps", p_col="p_wald"),
                display=DisplayConfig(show_recombination=False, snp_labels=False),
                liftover=LiftoverConfig(lifter=lifter),
            )

    def test_validates_the_config_before_lifting(self, plotter, source_gwas_df):
        """A lead lost to liftover is not announced as auto-detected and then refused."""
        lifter = InMemoryLifter({("chr1", 1_999): 11_999, ("chr1", 2_999): 12_999})

        with warnings.catch_warnings():
            warnings.simplefilter("error")
            with pytest.raises(ValidationError, match="needs a lead"):
                plotter.plot(
                    source_gwas_df,
                    chrom=1,
                    start=500,
                    end=3_500,
                    columns=ColumnConfig(pos_col="ps", p_col="p_wald"),
                    ld=LDConfig(lead_pos=1_000, ld_reference_file="/no/such/ref"),
                    display=DisplayConfig(show_recombination=False),
                    liftover=LiftoverConfig(lifter=lifter),
                )

    def test_an_invalid_region_is_refused_before_lifting(self, plotter, source_gwas_df):
        with pytest.raises(ValueError, match="must be < end"):
            plotter.plot(
                source_gwas_df,
                chrom=1,
                start=3_500,
                end=500,
                columns=ColumnConfig(pos_col="ps", p_col="p_wald"),
                liftover=LiftoverConfig(lifter=InMemoryLifter({("chr1", 0): 0})),
            )

    def test_plot_stacked_lifts_every_panel_into_one_window(
        self, plotter, source_gwas_df
    ):
        second = source_gwas_df.assign(ps=[1_000, 2_000, 4_000])
        lifter = InMemoryLifter(
            {
                ("chr1", 999): 10_999,
                ("chr1", 1_999): 11_999,
                ("chr1", 2_999): 12_999,
                ("chr1", 3_999): 13_999,
            }
        )

        fig = plotter.plot_stacked(
            [source_gwas_df, second],
            chrom=1,
            start=500,
            end=4_500,
            columns=ColumnConfig(pos_col="ps", p_col="p_wald"),
            lead_positions=[1_000, 4_000],
            display=DisplayConfig(show_recombination=False, snp_labels=False),
            liftover=LiftoverConfig(lifter=lifter),
        )

        top, bottom = fig.axes[0], fig.axes[1]
        top_xs = {x for coll in top.collections for x, _ in coll.get_offsets()}
        bottom_xs = {x for coll in bottom.collections for x, _ in coll.get_offsets()}
        assert top_xs == {11_000, 12_000, 13_000}
        assert bottom_xs == {11_000, 12_000, 14_000}
        assert top.get_xlim() == (10_500, 14_500)

    def test_recombination_uses_the_same_lifter_when_asked(
        self, plotter, source_gwas_df, tmp_path
    ):
        """lift_recombination sends the plotter's maps through the caller's chain."""
        maps = tmp_path / "maps"
        maps.mkdir()
        (maps / "chr1_recomb.tsv").write_text(
            "chr\tpos\trate\tcM\n1\t1000\t0.5\t0.1\n1\t2000\t0.7\t0.2\n"
        )
        lifting = LocusZoomPlotter(
            species="canine",
            genome_build="canfam4",
            recomb_data_dir=str(maps),
        )

        def overlay_positions(lift_recombination):
            fig = lifting.plot(
                source_gwas_df,
                chrom=1,
                start=500,
                end=3_500,
                columns=ColumnConfig(pos_col="ps", p_col="p_wald"),
                display=DisplayConfig(snp_labels=False),
                liftover=LiftoverConfig(
                    lifter=self.LIFTER, lift_recombination=lift_recombination
                ),
            )
            return [x for ax in fig.axes[1:] for x in ax.get_lines()[0].get_xdata()]

        assert overlay_positions(False) == []
        assert overlay_positions(True) == [11_000, 12_000]


class TestRecombinationLifter:
    def test_caller_lifter_replaces_the_registered_chain(self, tmp_path):
        (tmp_path / "chr1_recomb.tsv").write_text(
            "chr\tpos\trate\tcM\n1\t1000\t0.5\t0.1\n1\t2000\t0.7\t0.2\n"
        )
        lifter = InMemoryLifter({("chr1", 999): 10_999, ("chr1", 1_999): 11_999})

        result = get_recombination_rate_for_region(
            chrom=1,
            start=10_000,
            end=12_500,
            species="canine",
            data_dir=str(tmp_path),
            genome_build="canfam4",
            lifter=lifter,
        )

        assert result["pos"].tolist() == [11_000, 12_000]
