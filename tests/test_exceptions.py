"""Tests for the exception hierarchy and where the specialised errors are raised."""

import pandas as pd
import pytest

import pylocuszoom
from pylocuszoom.exceptions import (
    DataDownloadError,
    EmptyLDOutputError,
    EnsemblAPIError,
    EQTLValidationError,
    FinemappingValidationError,
    ForestValidationError,
    LoaderValidationError,
    OptionalDependencyMissing,
    PheWASValidationError,
    PlinkError,
    PyLocusZoomError,
    ReferenceAPIError,
    UCSCAPIError,
    ValidationError,
)

HIERARCHY = [
    (PyLocusZoomError, (Exception,)),
    (ValidationError, (PyLocusZoomError, ValueError)),
    (EQTLValidationError, (ValidationError, PyLocusZoomError, ValueError)),
    (FinemappingValidationError, (ValidationError, PyLocusZoomError, ValueError)),
    (LoaderValidationError, (ValidationError, PyLocusZoomError, ValueError)),
    (PheWASValidationError, (ValidationError, PyLocusZoomError, ValueError)),
    (ForestValidationError, (ValidationError, PyLocusZoomError, ValueError)),
    (PlinkError, (PyLocusZoomError, RuntimeError)),
    (EmptyLDOutputError, (PlinkError, PyLocusZoomError, RuntimeError)),
    (OptionalDependencyMissing, (PyLocusZoomError, ImportError)),
    (DataDownloadError, (PyLocusZoomError, RuntimeError)),
    (ReferenceAPIError, (DataDownloadError, PyLocusZoomError, RuntimeError)),
    (EnsemblAPIError, (ReferenceAPIError, DataDownloadError, PyLocusZoomError)),
    (UCSCAPIError, (ReferenceAPIError, DataDownloadError, PyLocusZoomError)),
]
"""Each exception and every class a caller may catch it as."""

EXCEPTIONS = [cls for cls, _ in HIERARCHY]


class TestExceptionHierarchy:
    """Test that all exceptions have correct inheritance."""

    @pytest.mark.parametrize(
        ("cls", "bases"), HIERARCHY, ids=[cls.__name__ for cls in EXCEPTIONS]
    )
    def test_exception_is_catchable_as_each_base(self, cls, bases):
        for base in bases:
            with pytest.raises(base):
                raise cls("test")

    def test_reference_api_error_is_a_download_error(self):
        """A service failure is a download failure, not an input validation one."""
        assert not issubclass(ReferenceAPIError, ValidationError)
        assert not issubclass(ReferenceAPIError, ValueError)

    @pytest.mark.parametrize("cls", EXCEPTIONS, ids=lambda cls: cls.__name__)
    def test_exception_message_preserved(self, cls):
        """The message a raise site passes is what str() shows."""
        assert str(cls("Column 'pos' is missing")) == "Column 'pos' is missing"

    @pytest.mark.parametrize("cls", EXCEPTIONS, ids=lambda cls: cls.__name__)
    def test_exceptions_importable_from_package(self, cls):
        """Every exception is importable from the top-level package."""
        assert getattr(pylocuszoom, cls.__name__) is cls


class TestExceptionChaining:
    """Test that exception chaining works correctly."""

    def test_raise_from_preserves_cause(self):
        """Exception chaining with 'raise X from Y' preserves __cause__."""
        original = ValueError("original error")
        with pytest.raises(ValidationError) as exc_info:
            try:
                raise original
            except ValueError as e:
                raise ValidationError("wrapped error") from e

        assert exc_info.value.__cause__ is original


class TestSpecializedExceptionsInUse:
    """The PheWAS and forest validators raise their own subclasses."""

    def test_phewas_validation_raises_phewas_error(self):
        """validate_phewas_df raises PheWASValidationError, not generic ValidationError."""
        from pylocuszoom.schemas import validate_phewas_df

        with pytest.raises(PheWASValidationError):
            validate_phewas_df(pd.DataFrame({"wrong_col": [1]}))

    def test_forest_validation_raises_forest_error(self):
        """validate_forest_df raises ForestValidationError, not generic ValidationError."""
        from pylocuszoom.schemas import validate_forest_df

        with pytest.raises(ForestValidationError):
            validate_forest_df(pd.DataFrame({"wrong_col": [1]}))
