# tests/test_ensembl.py
"""Tests for Ensembl REST API integration."""

import pytest


def test_get_ensembl_species_name_canine():
    """Test species name mapping for canine."""
    from pylocuszoom.ensembl import get_ensembl_species_name

    assert get_ensembl_species_name("canine") == "canis_lupus_familiaris"
    assert get_ensembl_species_name("dog") == "canis_lupus_familiaris"


def test_get_ensembl_species_name_human():
    """Test species name mapping for human."""
    from pylocuszoom.ensembl import get_ensembl_species_name

    assert get_ensembl_species_name("human") == "homo_sapiens"
    assert get_ensembl_species_name("homo_sapiens") == "homo_sapiens"


def test_get_ensembl_species_name_unknown():
    """Test unknown species returns input unchanged."""
    from pylocuszoom.ensembl import get_ensembl_species_name

    assert get_ensembl_species_name("my_custom_species") == "my_custom_species"
