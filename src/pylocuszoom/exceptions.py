"""Exception hierarchy for pyLocusZoom.

All pyLocusZoom exceptions inherit from PyLocusZoomError, enabling users to
catch all library errors with `except PyLocusZoomError`.
"""


class PyLocusZoomError(Exception):
    """Base exception for all pyLocusZoom errors."""


class ValidationError(PyLocusZoomError, ValueError):
    """Raised when input validation fails. Inherits ValueError for backward compat."""


class EQTLValidationError(ValidationError):
    """Raised when eQTL DataFrame validation fails."""


class FinemappingValidationError(ValidationError):
    """Raised when fine-mapping DataFrame validation fails."""


class LoaderValidationError(ValidationError):
    """Raised when loaded data fails validation."""


class PheWASValidationError(ValidationError):
    """Raised when PheWAS DataFrame validation fails."""


class ForestValidationError(ValidationError):
    """Raised when forest plot DataFrame validation fails."""


class PlinkError(PyLocusZoomError, RuntimeError):
    """Raised when PLINK subprocess fails."""


class EmptyLDOutputError(PlinkError):
    """Raised when PLINK succeeds but produces no LD pairs."""


class OptionalDependencyMissing(PyLocusZoomError, ImportError):
    """Raised when a feature needs an optional extra that is not installed."""


class DataDownloadError(PyLocusZoomError, RuntimeError):
    """Raised when data download operations fail."""


class ReferenceAPIError(DataDownloadError):
    """Raised when a reference-annotation API is unreachable or errors.

    Distinguishes a service failure from a rejected request, so callers can
    tell "the region has no genes" from "we could not ask" without knowing
    which source answered.
    """


class EnsemblAPIError(ReferenceAPIError):
    """Raised when the Ensembl REST API is unreachable or returns an error."""


class UCSCAPIError(ReferenceAPIError):
    """Raised when the UCSC REST API is unreachable or returns an error."""
