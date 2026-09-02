# tests/reference_mocks.py
"""Mock responses shared by the gene-source tests (Ensembl, UCSC, routing)."""

from unittest.mock import Mock


def ok_response(payload):
    """A mocked 200 response whose .json() returns the given payload."""
    response = Mock()
    response.ok = True
    response.status_code = 200
    response.json.return_value = payload
    return response


def ros_cfam_gene_payload():
    """ATP9B as Ensembl serves it for dog: ROS_Cfam_1.0 coordinates.

    CanFam3.1 puts the gene at 1,136,862; Ensembl answers a CanFam3.1
    request with these coordinates and an HTTP 200, which is the mismatch
    the assembly warning and the UCSC routing exist to catch.
    """
    return [
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


def gene_transcript_exon_payload():
    """One gene, one transcript and three exons, as Ensembl serves them.

    The transcript's ``Parent`` is the gene id and each exon's ``Parent`` is
    the transcript id, which is the join that gives an exon its gene symbol.
    """
    return [
        {
            "feature_type": "gene",
            "seq_region_name": "1",
            "start": 1_020_000,
            "end": 1_120_000,
            "external_name": "NFATC1",
            "strand": -1,
            "id": "ENSG01",
            "biotype": "protein_coding",
            "assembly_name": "GRCh38",
        },
        {
            "feature_type": "transcript",
            "seq_region_name": "1",
            "start": 1_020_000,
            "end": 1_120_000,
            "id": "ENST01",
            "Parent": "ENSG01",
            "assembly_name": "GRCh38",
        },
    ] + [
        {
            "feature_type": "exon",
            "seq_region_name": "1",
            "start": exon_start,
            "end": exon_start + 2_000,
            "id": f"ENSE0{index}",
            "Parent": "ENST01",
            "assembly_name": "GRCh38",
        }
        for index, exon_start in enumerate((1_020_000, 1_060_000, 1_110_000), start=1)
    ]


def refseq_payload():
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
