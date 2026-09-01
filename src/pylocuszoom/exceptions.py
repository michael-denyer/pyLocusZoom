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


class EnsemblAPIError(ValidationError):
    """Raised when the Ensembl REST API is unreachable or returns an error.

    Distinguishes a service failure from a rejected request, so callers can
    tell "the region has no genes" from "we could not ask". Subclasses
    ValidationError for backward compat with ``raise_on_error=True`` callers.
    """


class UCSCAPIError(ValidationError):
    """Raised when the UCSC REST API is unreachable or returns an error.

    The UCSC counterpart of EnsemblAPIError, with the same contract: a service
    failure stays distinguishable from a region that genuinely has no genes.
    """


class BackendError(PyLocusZoomError):
    """Raised when backend operations fail."""


class PlinkError(PyLocusZoomError, RuntimeError):
    """Raised when PLINK subprocess fails."""


class EmptyLDOutputError(PlinkError):
    """Raised when PLINK succeeds but produces no LD pairs."""


class DataDownloadError(PyLocusZoomError, RuntimeError):
    """Raised when data download operations fail."""
