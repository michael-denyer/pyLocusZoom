"""Tests for exception hierarchy."""

import pytest

from pylocuszoom.exceptions import (
    BackendError,
    DataDownloadError,
    EQTLValidationError,
    FinemappingValidationError,
    LoaderValidationError,
    PyLocusZoomError,
    ValidationError,
)


class TestExceptionHierarchy:
    """Test that all exceptions have correct inheritance."""

    def test_pylocuszoom_error_is_base_exception(self):
        """PyLocusZoomError inherits from Exception."""
        assert issubclass(PyLocusZoomError, Exception)

    def test_validation_error_inherits_from_base_and_value_error(self):
        """ValidationError inherits from both PyLocusZoomError and ValueError."""
        assert issubclass(ValidationError, PyLocusZoomError)
        assert issubclass(ValidationError, ValueError)

    def test_eqtl_validation_error_inherits_from_validation_error(self):
        """EQTLValidationError inherits from ValidationError."""
        assert issubclass(EQTLValidationError, ValidationError)
        assert issubclass(EQTLValidationError, PyLocusZoomError)
        assert issubclass(EQTLValidationError, ValueError)

    def test_finemapping_validation_error_inherits_from_validation_error(self):
        """FinemappingValidationError inherits from ValidationError."""
        assert issubclass(FinemappingValidationError, ValidationError)
        assert issubclass(FinemappingValidationError, PyLocusZoomError)
        assert issubclass(FinemappingValidationError, ValueError)

    def test_loader_validation_error_inherits_from_validation_error(self):
        """LoaderValidationError inherits from ValidationError."""
        assert issubclass(LoaderValidationError, ValidationError)
        assert issubclass(LoaderValidationError, PyLocusZoomError)
        assert issubclass(LoaderValidationError, ValueError)

    def test_backend_error_inherits_from_base(self):
        """BackendError inherits from PyLocusZoomError."""
        assert issubclass(BackendError, PyLocusZoomError)
        assert not issubclass(BackendError, ValueError)

    def test_data_download_error_inherits_from_base_and_runtime_error(self):
        """DataDownloadError inherits from both PyLocusZoomError and RuntimeError."""
        assert issubclass(DataDownloadError, PyLocusZoomError)
        assert issubclass(DataDownloadError, RuntimeError)


class TestExceptionInstantiation:
    """Test that all exceptions can be instantiated with messages."""

    def test_pylocuszoom_error_with_message(self):
        """PyLocusZoomError can be instantiated with message."""
        err = PyLocusZoomError("test message")
        assert str(err) == "test message"

    def test_validation_error_with_message(self):
        """ValidationError can be instantiated with message."""
        err = ValidationError("validation failed")
        assert str(err) == "validation failed"

    def test_eqtl_validation_error_with_message(self):
        """EQTLValidationError can be instantiated with message."""
        err = EQTLValidationError("eQTL validation failed")
        assert str(err) == "eQTL validation failed"

    def test_finemapping_validation_error_with_message(self):
        """FinemappingValidationError can be instantiated with message."""
        err = FinemappingValidationError("finemapping validation failed")
        assert str(err) == "finemapping validation failed"

    def test_loader_validation_error_with_message(self):
        """LoaderValidationError can be instantiated with message."""
        err = LoaderValidationError("loader validation failed")
        assert str(err) == "loader validation failed"

    def test_backend_error_with_message(self):
        """BackendError can be instantiated with message."""
        err = BackendError("backend failed")
        assert str(err) == "backend failed"

    def test_data_download_error_with_message(self):
        """DataDownloadError can be instantiated with message."""
        err = DataDownloadError("download failed")
        assert str(err) == "download failed"


class TestExceptionChaining:
    """Test that exception chaining works correctly."""

    def test_raise_from_preserves_cause(self):
        """Exception chaining with 'raise X from Y' preserves __cause__."""
        original = ValueError("original error")
        try:
            try:
                raise original
            except ValueError as e:
                raise ValidationError("wrapped error") from e
        except ValidationError as err:
            assert err.__cause__ is original
            assert str(err.__cause__) == "original error"

    def test_eqtl_error_wraps_validation_error(self):
        """EQTLValidationError can wrap ValidationError."""
        original = ValidationError("original validation error")
        try:
            try:
                raise original
            except ValidationError as e:
                raise EQTLValidationError("eQTL error") from e
        except EQTLValidationError as err:
            assert err.__cause__ is original


class TestCatchingExceptions:
    """Test that exceptions can be caught at various levels."""

    def test_catch_validation_error_as_value_error(self):
        """ValidationError can be caught as ValueError for backward compat."""
        with pytest.raises(ValueError):
            raise ValidationError("test")

    def test_catch_eqtl_error_as_validation_error(self):
        """EQTLValidationError can be caught as ValidationError."""
        with pytest.raises(ValidationError):
            raise EQTLValidationError("test")

    def test_catch_eqtl_error_as_value_error(self):
        """EQTLValidationError can be caught as ValueError."""
        with pytest.raises(ValueError):
            raise EQTLValidationError("test")

    def test_catch_eqtl_error_as_base(self):
        """EQTLValidationError can be caught as PyLocusZoomError."""
        with pytest.raises(PyLocusZoomError):
            raise EQTLValidationError("test")

    def test_catch_all_library_errors_as_base(self):
        """All library exceptions can be caught as PyLocusZoomError."""
        exceptions = [
            ValidationError("v"),
            EQTLValidationError("e"),
            FinemappingValidationError("f"),
            LoaderValidationError("l"),
            BackendError("b"),
            DataDownloadError("d"),
        ]
        for exc in exceptions:
            with pytest.raises(PyLocusZoomError):
                raise exc
