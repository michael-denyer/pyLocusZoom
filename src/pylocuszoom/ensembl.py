# src/pylocuszoom/ensembl.py
"""Ensembl REST API integration for reference data fetching.

Provides functions to fetch gene and exon annotations from the Ensembl REST API
(https://rest.ensembl.org) for any species.

Note: Recombination rates are NOT available from Ensembl for most species.
Use species-specific recombination maps instead (see recombination.py).
"""

from .logging import logger

# Ensembl API limits regions to 5Mb
ENSEMBL_MAX_REGION_SIZE = 5_000_000

# Species name aliases -> Ensembl species names
SPECIES_ALIASES: dict[str, str] = {
    # Canine
    "canine": "canis_lupus_familiaris",
    "dog": "canis_lupus_familiaris",
    "canis_familiaris": "canis_lupus_familiaris",
    # Feline
    "feline": "felis_catus",
    "cat": "felis_catus",
    # Human
    "human": "homo_sapiens",
    # Mouse
    "mouse": "mus_musculus",
    # Rat
    "rat": "rattus_norvegicus",
}


def get_ensembl_species_name(species: str) -> str:
    """Convert species alias to Ensembl species name.

    Args:
        species: Species name or alias (e.g., "canine", "dog", "human").

    Returns:
        Ensembl-compatible species name (e.g., "canis_lupus_familiaris").
    """
    return SPECIES_ALIASES.get(species.lower(), species.lower())
