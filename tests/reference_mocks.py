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
