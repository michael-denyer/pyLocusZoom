"""Exception hierarchy for pyLocusZoom.

All pyLocusZoom exceptions inherit from PyLocusZoomError, enabling users to
catch all library errors with `except PyLocusZoomError`.
"""


# [5a:PyLocusZoomError] Exception hierarchy root — see docs/CODEMAP.md
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


class BackendError(PyLocusZoomError):
    """Raised when backend operations fail."""


class PlinkError(PyLocusZoomError, RuntimeError):
    """Raised when PLINK subprocess fails."""


class DataDownloadError(PyLocusZoomError, RuntimeError):
    """Raised when data download operations fail."""
