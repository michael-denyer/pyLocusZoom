# tests/test_ucsc.py
"""Tests for the UCSC gene source and the build-to-source router.

Ensembl retired CanFam3.1, CanFam4 and FelCat9 and its archive REST hosts
redirect to a help page, so these builds have no Ensembl source at any URL.
UCSC hosts all three, which is the only reason this client exists.
"""

from unittest.mock import Mock, patch

import pandas as pd
import pytest


def _ok_response(payload):
    response = Mock()
    response.ok = True
    response.status_code = 200
    response.json.return_value = payload
    return response


def _refseq_payload():
    """Two NFATC1 transcripts and one non-coding gene, as UCSC returns them."""
    return {
        "ncbiRefSeq": [
            {
                "name": "XM_005615289.3",
                "name2": "NFATC1",
                "chrom": "chr1",
                "strand": "-",
                "txStart": 1027020,
                "txEnd": 1121187,
                "exonStarts": "1027020,1060176,",
                "exonEnds": "1027500,1060900,",
            },
            {
                "name": "XM_022408507.1",
                "name2": "NFATC1",
                "chrom": "chr1",
                "strand": "-",
                "txStart": 1060176,
                "txEnd": 1121362,
                "exonStarts": "1060176,",
                "exonEnds": "1060900,",
            },
            {
                "name": "XR_002615165.1",
                "name2": "LOC111090558",
                "chrom": "chr1",
                "strand": "+",
                "txStart": 1003647,
                "txEnd": 1004169,
                "exonStarts": "1003647,",
                "exonEnds": "1004169,",
            },
        ]
    }


class TestUCSCGeneFetch:
    def test_transcripts_collapse_to_one_gene(self):
        """Two transcripts of one symbol become one row spanning both."""
        from pylocuszoom.ucsc import fetch_genes_from_ucsc

        with patch(
            "pylocuszoom.ucsc.requests.get",
            return_value=_ok_response(_refseq_payload()),
        ):
            genes = fetch_genes_from_ucsc("canFam3", "1", 1_000_000, 1_200_000)

        assert genes["gene_name"].tolist() == ["NFATC1"]
        row = genes.iloc[0]
        assert row["start"] == 1027021, "UCSC txStart is 0-based, output is 1-based"
        assert row["end"] == 1121362, "widest transcript end wins"
        assert row["strand"] == "-"

    def test_biotype_filter_drops_non_coding(self):
        """The default protein_coding filter drops XR_/NR_ only genes."""
        from pylocuszoom.ucsc import fetch_genes_from_ucsc

        with patch(
            "pylocuszoom.ucsc.requests.get",
            return_value=_ok_response(_refseq_payload()),
        ):
            all_genes = fetch_genes_from_ucsc(
                "canFam3", "1", 1_000_000, 1_200_000, biotype=""
            )

        assert set(all_genes["gene_name"]) == {"NFATC1", "LOC111090558"}
        non_coding = all_genes[all_genes.gene_name == "LOC111090558"]
        assert non_coding["biotype"].iloc[0] == "non_coding"

    def test_rows_record_the_assembly(self):
        """Every row names the assembly, matching the Ensembl frame's shape."""
        from pylocuszoom.ucsc import fetch_genes_from_ucsc

        with patch(
            "pylocuszoom.ucsc.requests.get",
            return_value=_ok_response(_refseq_payload()),
        ):
            genes = fetch_genes_from_ucsc("canFam3", "1", 1_000_000, 1_200_000)

        assert genes["assembly"].tolist() == ["CanFam3.1"]

    def test_exons_carry_their_gene_name(self):
        """Exons name their gene, so the gene track can match them to a row."""
        from pylocuszoom.ucsc import fetch_exons_from_ucsc

        with patch(
            "pylocuszoom.ucsc.requests.get",
            return_value=_ok_response(_refseq_payload()),
        ):
            exons = fetch_exons_from_ucsc("canFam3", "1", 1_000_000, 1_200_000)

        assert set(exons["gene_name"]) == {"NFATC1", "LOC111090558"}
        assert exons["start"].min() == 1003648, "exonStarts are 0-based too"
        assert exons["assembly"].unique().tolist() == ["CanFam3.1"]

    def test_api_error_returns_empty_frame_with_columns(self):
        """An outage yields an empty frame that still round-trips through CSV."""
        import requests

        from pylocuszoom.ucsc import fetch_genes_from_ucsc

        with (
            patch("pylocuszoom.ucsc.time.sleep"),
            patch(
                "pylocuszoom.ucsc.requests.get",
                side_effect=requests.exceptions.ConnectionError("network down"),
            ),
        ):
            genes = fetch_genes_from_ucsc("canFam3", "1", 1_000_000, 1_200_000)

        assert genes.empty
        assert "gene_name" in genes.columns
        assert "assembly" in genes.columns

    def test_api_error_raises_when_asked(self):
        """raise_on_error keeps a service failure distinct from an empty region."""
        import requests

        from pylocuszoom.exceptions import UCSCAPIError
        from pylocuszoom.ucsc import fetch_genes_from_ucsc

        with (
            patch("pylocuszoom.ucsc.time.sleep"),
            patch(
                "pylocuszoom.ucsc.requests.get",
                side_effect=requests.exceptions.ConnectionError("network down"),
            ),
            pytest.raises(UCSCAPIError),
        ):
            fetch_genes_from_ucsc(
                "canFam3", "1", 1_000_000, 1_200_000, raise_on_error=True
            )


class TestUCSCCaching:
    def test_second_call_is_served_from_cache(self, tmp_path):
        from pylocuszoom.ucsc import get_genes_for_region_ucsc

        with patch(
            "pylocuszoom.ucsc.requests.get",
            return_value=_ok_response(_refseq_payload()),
        ) as mock_get:
            first = get_genes_for_region_ucsc(
                "canFam3", "1", 1_000_000, 1_200_000, tmp_path
            )
            second = get_genes_for_region_ucsc(
                "canFam3", "1", 1_000_000, 1_200_000, tmp_path
            )

        assert mock_get.call_count == 1
        assert second["gene_name"].tolist() == first["gene_name"].tolist()

    def test_failed_fetch_is_not_cached(self, tmp_path):
        """An outage must not poison the region with a permanent empty result."""
        import requests

        from pylocuszoom.ucsc import get_genes_for_region_ucsc

        with (
            patch("pylocuszoom.ucsc.time.sleep"),
            patch(
                "pylocuszoom.ucsc.requests.get",
                side_effect=requests.exceptions.ConnectionError("network down"),
            ),
        ):
            during_outage = get_genes_for_region_ucsc(
                "canFam3", "1", 1_000_000, 1_200_000, tmp_path
            )
        assert during_outage.empty

        with patch(
            "pylocuszoom.ucsc.requests.get",
            return_value=_ok_response(_refseq_payload()),
        ):
            after = get_genes_for_region_ucsc(
                "canFam3", "1", 1_000_000, 1_200_000, tmp_path
            )
        assert after["gene_name"].tolist() == ["NFATC1"]

    def test_two_genomes_do_not_share_an_entry(self, tmp_path):
        from pylocuszoom.ucsc import get_genes_for_region_ucsc

        with patch(
            "pylocuszoom.ucsc.requests.get",
            return_value=_ok_response(_refseq_payload()),
        ) as mock_get:
            get_genes_for_region_ucsc("canFam3", "1", 1_000_000, 1_200_000, tmp_path)
            get_genes_for_region_ucsc("canFam4", "1", 1_000_000, 1_200_000, tmp_path)

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
            "pylocuszoom.ucsc.requests.get",
            return_value=_ok_response(_refseq_payload()),
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

        ensembl_payload = [
            {
                "feature_type": "gene",
                "seq_region_name": "1",
                "start": 938796,
                "end": 1175952,
                "external_name": "ATP9B",
                "strand": -1,
                "id": "ENSCAFG00845000134",
                "biotype": "protein_coding",
                "assembly_name": "ROS_Cfam_1.0",
            }
        ]
        with patch(
            "pylocuszoom.ensembl.requests.get",
            return_value=_ok_response(ensembl_payload),
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
            return_value=pd.DataFrame(
                columns=["chr", "start", "end", "gene_name", "assembly"]
            ),
        ) as mock_fetch:
            plotter.plot(gwas, chrom=1, start=900_000, end=1_200_000)

        assert mock_fetch.call_args.kwargs["genome_build"] == plotter.genome_build
