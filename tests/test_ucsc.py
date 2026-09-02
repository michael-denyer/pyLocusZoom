# tests/test_ucsc.py
"""Tests for the UCSC gene source and the build-to-source router.

Ensembl retired CanFam3.1, CanFam4 and FelCat9 and its archive REST hosts
redirect to a help page, so these builds have no Ensembl source at any URL.
UCSC hosts all three, which is the only reason this client exists.
"""

from unittest.mock import patch

import pandas as pd
import pytest

from tests.reference_mocks import ok_response, refseq_payload, ros_cfam_gene_payload


class TestUCSCGeneFetch:
    def test_transcripts_collapse_to_one_gene(self):
        """Two transcripts of one symbol become one row spanning both."""
        from pylocuszoom.ucsc import fetch_track_frames

        with patch(
            "pylocuszoom._http.requests.get",
            return_value=ok_response(refseq_payload()),
        ):
            genes, _ = fetch_track_frames("canFam3", "1", 1_000_000, 1_200_000)

        assert genes["gene_name"].tolist() == ["NFATC1"]
        row = genes.iloc[0]
        assert row["start"] == 1027021, "UCSC txStart is 0-based, output is 1-based"
        assert row["end"] == 1121362, "widest transcript end wins"
        assert row["strand"] == "-"

    def test_biotype_filter_drops_non_coding(self):
        """The default protein_coding filter drops XR_/NR_ only genes."""
        from pylocuszoom.ucsc import fetch_track_frames

        with patch(
            "pylocuszoom._http.requests.get",
            return_value=ok_response(refseq_payload()),
        ):
            all_genes, _ = fetch_track_frames(
                "canFam3", "1", 1_000_000, 1_200_000, biotype=""
            )

        assert set(all_genes["gene_name"]) == {"NFATC1", "LOC111090558"}
        non_coding = all_genes[all_genes.gene_name == "LOC111090558"]
        assert non_coding["biotype"].iloc[0] == "non_coding"

    def test_rows_record_the_assembly(self):
        """Every row names the assembly, matching the Ensembl frame's shape."""
        from pylocuszoom.ucsc import fetch_track_frames

        with patch(
            "pylocuszoom._http.requests.get",
            return_value=ok_response(refseq_payload()),
        ):
            genes, _ = fetch_track_frames("canFam3", "1", 1_000_000, 1_200_000)

        assert genes["assembly"].tolist() == ["CanFam3.1"]

    def test_exons_carry_their_gene_name(self):
        """Exons name their gene, so the gene track can match them to a row."""
        from pylocuszoom.ucsc import fetch_track_frames

        with patch(
            "pylocuszoom._http.requests.get",
            return_value=ok_response(refseq_payload()),
        ):
            _, exons = fetch_track_frames("canFam3", "1", 1_000_000, 1_200_000)

        assert set(exons["gene_name"]) == {"NFATC1", "LOC111090558"}
        assert exons["start"].min() == 1003648, "exonStarts are 0-based too"
        assert exons["assembly"].unique().tolist() == ["CanFam3.1"]

    def test_empty_region_keeps_the_columns(self):
        """A gene-free region yields frames that still round-trip through CSV."""
        from pylocuszoom.ucsc import fetch_track_frames

        with patch(
            "pylocuszoom._http.requests.get",
            return_value=ok_response({"ncbiRefSeq": []}),
        ):
            genes, exons = fetch_track_frames("canFam3", "1", 1_000_000, 1_200_000)

        assert genes.empty and exons.empty
        assert "gene_name" in genes.columns
        assert "assembly" in genes.columns
        assert "exon_id" in exons.columns

    def test_api_error_raises(self):
        """A service failure stays distinct from an empty region."""
        import requests

        from pylocuszoom.exceptions import UCSCAPIError
        from pylocuszoom.ucsc import fetch_track_frames

        with (
            patch("pylocuszoom._http.time.sleep"),
            patch(
                "pylocuszoom._http.requests.get",
                side_effect=requests.exceptions.ConnectionError("network down"),
            ),
            pytest.raises(UCSCAPIError),
        ):
            fetch_track_frames("canFam3", "1", 1_000_000, 1_200_000)


class TestUCSCCaching:
    def test_second_call_is_served_from_cache(self, tmp_path):
        from pylocuszoom.reference_genes import get_genes_for_build
        from pylocuszoom.ucsc import ucsc_source

        with patch(
            "pylocuszoom._http.requests.get",
            return_value=ok_response(refseq_payload()),
        ) as mock_get:
            first = get_genes_for_build(
                "",
                "1",
                1_000_000,
                1_200_000,
                cache_dir=tmp_path,
                source=ucsc_source("canFam3"),
            )
            second = get_genes_for_build(
                "",
                "1",
                1_000_000,
                1_200_000,
                cache_dir=tmp_path,
                source=ucsc_source("canFam3"),
            )

        assert mock_get.call_count == 1
        assert second["gene_name"].tolist() == first["gene_name"].tolist()

    def test_failed_fetch_is_not_cached(self, tmp_path):
        """An outage must not poison the region with a permanent empty result."""
        import requests

        from pylocuszoom.exceptions import UCSCAPIError
        from pylocuszoom.reference_genes import get_genes_for_build
        from pylocuszoom.ucsc import ucsc_source

        with (
            patch("pylocuszoom._http.time.sleep"),
            patch(
                "pylocuszoom._http.requests.get",
                side_effect=requests.exceptions.ConnectionError("network down"),
            ),
            pytest.raises(UCSCAPIError),
        ):
            get_genes_for_build(
                "",
                "1",
                1_000_000,
                1_200_000,
                cache_dir=tmp_path,
                source=ucsc_source("canFam3"),
            )

        with patch(
            "pylocuszoom._http.requests.get",
            return_value=ok_response(refseq_payload()),
        ):
            after = get_genes_for_build(
                "",
                "1",
                1_000_000,
                1_200_000,
                cache_dir=tmp_path,
                source=ucsc_source("canFam3"),
            )
        assert after["gene_name"].tolist() == ["NFATC1"]

    def test_exons_share_the_gene_request(self, tmp_path):
        """ncbiRefSeq rows carry exons, so one track request serves both."""
        from pylocuszoom.reference_genes import get_genes_for_build
        from pylocuszoom.ucsc import ucsc_source

        with patch(
            "pylocuszoom._http.requests.get",
            return_value=ok_response(refseq_payload()),
        ) as mock_get:
            genes, exons = get_genes_for_build(
                "",
                "1",
                1_000_000,
                1_200_000,
                cache_dir=tmp_path,
                include_exons=True,
                source=ucsc_source("canFam3"),
            )

        assert mock_get.call_count == 1
        assert genes["gene_name"].tolist() == ["NFATC1"]
        assert set(exons["gene_name"]) == {"NFATC1", "LOC111090558"}

    def test_two_genomes_do_not_share_an_entry(self, tmp_path):
        from pylocuszoom.reference_genes import get_genes_for_build
        from pylocuszoom.ucsc import ucsc_source

        with patch(
            "pylocuszoom._http.requests.get",
            return_value=ok_response(refseq_payload()),
        ) as mock_get:
            for genome in ("canFam3", "canFam4"):
                get_genes_for_build(
                    "",
                    "1",
                    1_000_000,
                    1_200_000,
                    cache_dir=tmp_path,
                    source=ucsc_source(genome),
                )

        assert mock_get.call_count == 2


class TestBuildRouting:
    """The build decides the source, and that decision is the whole policy."""

    @pytest.mark.parametrize(
        "build,expected",
        [
            ("canfam3.1", "canFam3"),
            ("CanFam3", "canFam3"),
            ("canfam4", "canFam4"),
            ("UU_Cfam_GSD_1.0", "canFam4"),
            ("felCat9", "felCat9"),
            ("Felis_catus_9.0", "felCat9"),
            ("ROS_Cfam_1.0", None),
            ("GRCh38", None),
            (None, None),
        ],
    )
    def test_source_selection(self, build, expected):
        from pylocuszoom.reference_genes import ucsc_genome_for_build

        assert ucsc_genome_for_build(build) == expected

    def test_retired_build_routes_to_ucsc(self, tmp_path):
        """A CanFam3.1 caller gets CanFam3.1 coordinates, not ROS_Cfam_1.0."""
        from pylocuszoom.reference_genes import get_genes_for_build

        with patch(
            "pylocuszoom._http.requests.get",
            return_value=ok_response(refseq_payload()),
        ):
            genes = get_genes_for_build(
                "canine",
                "1",
                1_000_000,
                1_200_000,
                genome_build="canfam3.1",
                cache_dir=tmp_path,
            )

        assert genes["assembly"].tolist() == ["CanFam3.1"]
        assert genes.iloc[0]["start"] == 1027021

    def test_served_build_stays_on_ensembl(self, tmp_path):
        """A build Ensembl serves is not diverted to UCSC."""
        from pylocuszoom.reference_genes import get_genes_for_build

        with patch(
            "pylocuszoom._http.requests.get",
            return_value=ok_response(ros_cfam_gene_payload()),
        ):
            genes = get_genes_for_build(
                "canine",
                "1",
                900_000,
                1_200_000,
                genome_build="ROS_Cfam_1.0",
                cache_dir=tmp_path,
            )

        assert genes["assembly"].tolist() == ["ROS_Cfam_1.0"]

    def test_plotter_default_canine_build_reaches_ucsc(self):
        """The canine default is CanFam3.1, so auto_genes must land on UCSC."""
        from pylocuszoom import LocusZoomPlotter
        from pylocuszoom.reference_genes import ucsc_genome_for_build

        plotter = LocusZoomPlotter(species="canine", auto_genes=True, log_level=None)
        assert ucsc_genome_for_build(plotter.genome_build) == "canFam3"

    def test_plotter_fetches_in_its_own_build(self):
        """auto_genes routes through the build-aware entry point."""
        from pylocuszoom import LocusZoomPlotter

        plotter = LocusZoomPlotter(species="canine", auto_genes=True, log_level=None)
        gwas = pd.DataFrame(
            {"ps": [1_000_000, 1_050_000], "p_wald": [0.5, 1e-6], "rs": ["a", "b"]}
        )

        with patch(
            "pylocuszoom.plotter.get_genes_for_build",
            return_value=(
                pd.DataFrame(columns=["chr", "start", "end", "gene_name", "assembly"]),
                pd.DataFrame(),
            ),
        ) as mock_fetch:
            plotter.plot(gwas, chrom=1, start=900_000, end=1_200_000)

        assert mock_fetch.call_args.kwargs["genome_build"] == plotter.genome_build
