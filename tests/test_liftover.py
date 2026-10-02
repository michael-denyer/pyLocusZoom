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
from hypothesis import assume, given
from hypothesis import strategies as st
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
from pylocuszoom._liftover import InMemoryLifter, chain_lifter, lift_window
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
    """Registered chains download once into the cache, and only if they parse."""

    CANFAM3, CANFAM4 = GENOME_BUILDS["canfam3"], GENOME_BUILDS["canfam4"]

    @pytest.fixture
    def downloads(self, monkeypatch):
        """Record every chain download, writing SPLIT_CHAIN to the destination."""
        urls = []

        def download(url, dest, desc=None):
            urls.append(url)
            content = SPLIT_CHAIN.encode()
            dest.write_bytes(gzip.compress(content) if url.endswith(".gz") else content)

        monkeypatch.setattr("pylocuszoom._liftover.stream_file", download)
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

        lifter = chain_lifter(source, target)

        assert downloads == [url]
        assert lifter.convert_coordinate("chr1", 999)[0][:2] == ("chr1", 1099)

    def test_an_unregistered_pair_is_a_validation_error(self, cache_home, downloads):
        with pytest.raises(ValidationError, match="No liftover chain"):
            chain_lifter(self.CANFAM4, self.CANFAM3)
        assert downloads == []

    def test_an_unreadable_download_is_reported_and_not_cached(
        self, cache_home, monkeypatch
    ):
        def download(url, dest, desc=None):
            dest.write_bytes(b"<html>502</html>")

        monkeypatch.setattr("pylocuszoom._liftover.stream_file", download)

        with pytest.raises(DataDownloadError, match="unreadable") as raised:
            chain_lifter(self.CANFAM3, self.CANFAM4)
        assert not (cache_home / "liftover" / "canFam3ToCanFam4.over.chain.gz").exists()
        # The staged file is gone by now, so the message names the URL once.
        message = str(raised.value)
        assert self.CANFAM3.chain_url(self.CANFAM4) in message
        assert ".part" not in message
        assert message.count("is unreadable") == 1

    def test_a_chain_download_is_staged_once(self, cache_home, monkeypatch):
        """The download is written straight into the sibling that gets parsed."""
        written = []

        def stream(url, partial, desc=None):
            written.append(sorted(p.name for p in partial.parent.iterdir()))
            partial.write_bytes(gzip.compress(SPLIT_CHAIN.encode()))

        monkeypatch.setattr("pylocuszoom._liftover.stream_file", stream)

        chain_lifter(self.CANFAM3, self.CANFAM4)

        assert len(written) == 1 and len(written[0]) == 1
        assert [p.name for p in (cache_home / "liftover").iterdir()] == [
            "canFam3ToCanFam4.over.chain.gz"
        ]

    CORRUPT = b"<html>captive portal</html>"

    @pytest.fixture
    def chain_path(self, cache_home):
        return cache_home / "liftover" / "canFam3ToCanFam4.over.chain.gz"

    @pytest.fixture
    def serve(self, monkeypatch):
        """Make successive chain downloads write the given bodies, or raise them."""

        def serve(*bodies):
            remaining = list(bodies)

            def download(url, dest, desc=None):
                body = remaining.pop(0)
                if isinstance(body, Exception):
                    raise body
                dest.write_bytes(body)

            monkeypatch.setattr("pylocuszoom._liftover.stream_file", download)

        return serve

    def test_a_corrupt_download_leaves_another_process_chain_in_place(
        self, chain_path, serve, monkeypatch
    ):
        """Another process publishes a good chain when this one fails to load."""
        from pylocuszoom import _liftover

        good = gzip.compress(SPLIT_CHAIN.encode())
        chain_path.parent.mkdir(parents=True)
        chain_path.write_bytes(self.CORRUPT)
        serve(self.CORRUPT)
        real_load = _liftover.load_chain

        def load(path):
            try:
                return real_load(path)
            except ValidationError:
                chain_path.write_bytes(good)
                raise

        monkeypatch.setattr(_liftover, "load_chain", load)

        with pytest.raises(DataDownloadError, match="unreadable"):
            chain_lifter(self.CANFAM3, self.CANFAM4)
        assert chain_path.read_bytes() == good

    def test_a_corrupt_download_is_never_cached_while_downloads_keep_failing(
        self, chain_path, serve
    ):
        serve(self.CORRUPT, DataDownloadError("network down"))

        for _ in range(2):
            with pytest.raises(DataDownloadError):
                chain_lifter(self.CANFAM3, self.CANFAM4)
            assert list(chain_path.parent.iterdir()) == []

    def test_a_held_lifter_is_served_after_the_cached_chain_is_removed(
        self, chain_path, serve
    ):
        serve(gzip.compress(SPLIT_CHAIN.encode()), DataDownloadError("network down"))
        first = chain_lifter(self.CANFAM3, self.CANFAM4)
        chain_path.unlink()

        assert chain_lifter(self.CANFAM3, self.CANFAM4) is first

    def test_a_corrupt_cached_chain_is_replaced_by_a_good_download(
        self, chain_path, serve
    ):
        good = gzip.compress(SPLIT_CHAIN.encode())
        chain_path.parent.mkdir(parents=True)
        chain_path.write_bytes(self.CORRUPT)
        serve(good)

        lifter = chain_lifter(self.CANFAM3, self.CANFAM4)

        assert lifter.convert_coordinate("chr1", 999)[0][:2] == ("chr1", 1099)
        assert chain_path.read_bytes() == good

    def test_a_corrupt_cached_chain_is_kept_when_the_download_is_corrupt(
        self, chain_path, serve
    ):
        """Removing it could remove a good chain another process just published."""
        chain_path.parent.mkdir(parents=True)
        chain_path.write_bytes(self.CORRUPT)
        serve(self.CORRUPT)

        with pytest.raises(DataDownloadError, match="unreadable"):
            chain_lifter(self.CANFAM3, self.CANFAM4)
        assert [p.name for p in chain_path.parent.iterdir()] == [chain_path.name]

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

    # 1-based 2500 is in the region but is no row of ``source_gwas_df``.
    STRAY_LEAD = 2_500

    def test_warns_when_lead_lifts_outside_the_window(self, plotter, source_gwas_df):
        lifter = InMemoryLifter({**self.LIFTER._mapping, ("chr1", 2_499): 99_999})
        with pytest.warns(
            UserWarning, match="Lead SNP at chr1:2500 lifted outside the window"
        ):
            fig = plotter.plot(
                source_gwas_df,
                chrom=1,
                start=500,
                end=3_500,
                columns=ColumnConfig(pos_col="ps", p_col="p_wald"),
                ld=LDConfig(lead_pos=self.STRAY_LEAD),
                display=DisplayConfig(show_recombination=False, snp_labels=False),
                liftover=LiftoverConfig(lifter=lifter),
            )

        assert fig.axes[0].get_xlim() == (10_500, 13_500)

    def test_warns_when_lead_lifts_onto_another_snp(self, plotter, source_gwas_df):
        """A lead that is no row must not turn the SNP it lands on into the lead."""
        lifter = InMemoryLifter({**self.LIFTER._mapping, ("chr1", 2_499): 11_999})
        with pytest.warns(
            UserWarning, match="Lead SNP at chr1:2500 is not a SNP of the data"
        ):
            plotter.plot(
                source_gwas_df,
                chrom=1,
                start=500,
                end=3_500,
                columns=ColumnConfig(pos_col="ps", p_col="p_wald"),
                ld=LDConfig(lead_pos=self.STRAY_LEAD),
                display=DisplayConfig(show_recombination=False, snp_labels=False),
                liftover=LiftoverConfig(lifter=lifter),
            )

    @pytest.mark.parametrize("lifted_lead", [99_999, 11_999])
    def test_an_unusable_lifted_lead_is_an_error_with_ld_reference_file(
        self, plotter, source_gwas_df, lifted_lead
    ):
        """Outside the window or on another SNP, PLINK is not run for a wrong lead."""
        lifter = InMemoryLifter({**self.LIFTER._mapping, ("chr1", 2_499): lifted_lead})

        with warnings.catch_warnings():
            warnings.simplefilter("error")
            with pytest.raises(ValidationError, match="chr1:2500 .* needs a lead"):
                plotter.plot(
                    source_gwas_df,
                    chrom=1,
                    start=500,
                    end=3_500,
                    columns=ColumnConfig(pos_col="ps", p_col="p_wald"),
                    ld=LDConfig(
                        lead_pos=self.STRAY_LEAD, ld_reference_file="/no/such/ref"
                    ),
                    display=DisplayConfig(show_recombination=False),
                    liftover=LiftoverConfig(lifter=lifter),
                )

    def test_plot_stacked_warns_when_a_lead_lifts_outside_the_window(
        self, plotter, source_gwas_df
    ):
        lifter = InMemoryLifter({**self.LIFTER._mapping, ("chr1", 2_499): 99_999})
        with pytest.warns(
            UserWarning, match="Lead SNP at chr1:2500 lifted outside the window"
        ):
            plotter.plot_stacked(
                [source_gwas_df, source_gwas_df],
                chrom=1,
                start=500,
                end=3_500,
                columns=ColumnConfig(pos_col="ps", p_col="p_wald"),
                lead_positions=[self.STRAY_LEAD, 1_000],
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


@pytest.mark.parametrize("content", [None, "chain bad header\n", "not a chain\n"])
def test_unreadable_caller_chain_is_a_validation_error(tmp_path, content):
    path = tmp_path / "invalid.chain"
    if content is not None:
        path.write_text(content)

    with pytest.raises(ValidationError, match="invalid.chain"):
        LiftoverConfig(chain_path=path).resolve()


def test_unwritable_chain_cache_is_a_download_error(cache_home):
    cache_home.mkdir(parents=True, exist_ok=True)
    (cache_home / "liftover").write_text("a file blocks creating the cache directory")

    with pytest.raises(DataDownloadError, match="liftover"):
        chain_lifter(GENOME_BUILDS["canfam3"], GENOME_BUILDS["canfam4"])


@st.composite
def _lift_cases(draw, identity=False):
    """A region, one to three frames around it, and a lifter for their SNPs."""
    start = draw(st.integers(min_value=1, max_value=1_000_000))
    end = start + draw(st.integers(min_value=1, max_value=1_000_000))
    frames, leads, mapping = [], [], {("chr1", -1): 0}
    for _ in range(draw(st.integers(min_value=1, max_value=3))):
        positions = draw(
            st.lists(
                st.integers(min_value=max(1, start - 1_000), max_value=end + 1_000),
                min_size=1,
                max_size=20,
            )
        )
        for pos in positions:
            target = (
                pos
                if identity
                else draw(st.none() | st.integers(min_value=1, max_value=5_000_000))
            )
            if target is not None:
                mapping[("chr1", pos - 1)] = target - 1
        frames.append(pd.DataFrame({"chr": 1, "pos": positions, "p_value": 0.5}))
        # A lead is any position of the region, a row of its frame or not.
        lead = draw(st.none() | st.integers(min_value=start, max_value=end))
        if lead is not None and not identity and ("chr1", lead - 1) not in mapping:
            # Off a row it lifts anywhere, or onto a position a SNP lifted to.
            target = draw(
                st.none()
                | st.integers(min_value=1, max_value=5_000_000)
                | st.sampled_from(sorted(mapping.values())).map(lambda hit: hit + 1)
            )
            if target is not None:
                mapping[("chr1", lead - 1)] = target - 1
        leads.append(lead)
    return start, end, frames, leads, InMemoryLifter(mapping)


def _lift(start, end, frames, leads, lifter):
    return lift_window(
        frames,
        chrom=1,
        start=start,
        end=end,
        chrom_col="chr",
        pos_col="pos",
        lead_positions=leads,
        lifter=lifter,
        build=None,
    )


class TestLiftWindowProperties:
    """Invariants of the window every regional liftover goes through."""

    @given(_lift_cases())
    def test_window_holds_every_lifted_snp_and_lead(self, case):
        try:
            window = _lift(*case)
        except ValidationError:
            return  # A frame with no lifted SNP in the region; pinned elsewhere.

        _, _, frames, leads, _ = case
        assert 1 <= window.start < window.end
        for frame in window.frames:
            assert window.start <= frame["pos"].min()
            assert frame["pos"].max() <= window.end
        for source, lifted, source_lead, lead in zip(
            frames, window.frames, leads, window.lead_positions
        ):
            if lead is None:
                continue
            assert window.start <= lead <= window.end
            # Only a lead that is a lifted row may share a lifted row's position.
            if source_lead not in source.loc[lifted.index, "pos"].tolist():
                assert lead not in lifted["pos"].tolist()

    @given(_lift_cases(identity=True))
    def test_an_identity_lift_keeps_the_requested_window(self, case):
        start, end, frames, _, _ = case
        assume(all(frame["pos"].between(start, end).any() for frame in frames))

        window = _lift(*case)

        assert (window.start, window.end) == (start, end)
        for before, after in zip(frames, window.frames):
            in_region = before[before["pos"].between(start, end)]
            assert after["pos"].tolist() == in_region["pos"].tolist()


class TestLiftWindowLeads:
    """A lead lifts on its own, so it is checked against the lifted rows."""

    @pytest.fixture
    def frame(self):
        return pd.DataFrame({"chr": 1, "pos": [100, 200], "p_value": [1e-9, 1e-3]})

    def test_a_lead_lifting_outside_the_window_is_dropped_with_a_note(self, frame):
        lifter = InMemoryLifter(
            {("chr1", 99): 99, ("chr1", 149): 4_999, ("chr1", 199): 199}
        )

        window = _lift(50, 250, [frame], [150], lifter)

        assert (window.start, window.end) == (50, 250)
        assert window.lead_positions == [None]
        assert window.notes == (
            "Lead SNP at chr1:150 lifted outside the window of the region's "
            "SNPs in the target build; the lead is auto-detected instead",
        )

    def test_a_lead_that_is_no_row_lifting_onto_a_snp_is_dropped_with_a_note(
        self, frame
    ):
        lifter = InMemoryLifter(
            {("chr1", 99): 99, ("chr1", 149): 299, ("chr1", 199): 299}
        )

        window = _lift(50, 250, [frame], [150], lifter)

        assert window.frames[0]["pos"].tolist() == [100, 300]
        assert window.lead_positions == [None]
        assert window.notes == (
            "Lead SNP at chr1:150 is not a SNP of the data and lifted onto "
            "another SNP's position in the target build; the lead is "
            "auto-detected instead",
        )

    def test_a_lead_that_is_a_row_keeps_its_lifted_position(self, frame):
        lifter = InMemoryLifter({("chr1", 99): 99, ("chr1", 199): 299})

        window = _lift(50, 250, [frame], [200], lifter)

        assert window.lead_positions == [300]
        assert window.notes == ()
