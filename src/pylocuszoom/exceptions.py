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


class BackendError(PyLocusZoomError):
    """Raised when backend operations fail."""


class DataDownloadError(PyLocusZoomError, RuntimeError):
    """Raised when data download operations fail."""
