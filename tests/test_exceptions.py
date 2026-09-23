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


_GWAS = pd.DataFrame(
    {"chr": ["1"] * 3, "pos": [100, 200, 300], "p_value": [0.1, 0.2, 0.3]}
)
_COLOC = pd.DataFrame(
    {"pos": [100, 200], "p_gwas": [0.1, 0.2], "p_eqtl": [0.1, 0.2], "rs": ["a", "b"]}
)


def _regional():
    return pylocuszoom.LocusZoomPlotter(species=None)


def _genomewide():
    return pylocuszoom.ManhattanPlotter(species="canine")


INPUT_ERRORS = [
    pytest.param(
        lambda: _regional().plot(_GWAS, chrom=1, start=500, end=100),
        "start",
        id="region start after end",
    ),
    pytest.param(
        lambda: _regional().plot_stacked([], chrom=1, start=1, end=1000),
        "GWAS DataFrame",
        id="plot_stacked empty list",
    ),
    pytest.param(
        lambda: _genomewide().plot_qq(_GWAS.assign(p_value=0.0)),
        "p-value",
        id="qq all-zero p",
    ),
    pytest.param(
        lambda: _genomewide().plot_manhattan(_GWAS.assign(p_value=float("nan"))),
        "p_value",
        id="manhattan all-NaN p",
    ),
    pytest.param(
        lambda: _genomewide().plot_manhattan(_GWAS, category_col="nope"),
        "nope",
        id="categorical missing column",
    ),
    pytest.param(
        lambda: _regional().plot(
            _GWAS.drop(columns="p_value"), chrom=1, start=1, end=1000
        ),
        "p_value",
        id="regional missing p column",
    ),
    pytest.param(
        lambda: pylocuszoom.LDConfig(ld_col="r2", ld_reference_file="x"),
        "ld_col",
        id="LDConfig both sources",
    ),
    pytest.param(
        lambda: pylocuszoom.GenomeWideStyle(palette=["not-a-colour"]),
        "palette",
        id="GenomeWideStyle bad palette",
    ),
    pytest.param(
        lambda: _genomewide().plot_qq(
            _GWAS, config=pylocuszoom.GenomeWideConfig(p_col="p")
        ),
        "'p'",
        id="qq missing p column",
    ),
    pytest.param(
        lambda: _genomewide().plot_manhattan_stacked([]),
        "GWAS DataFrame",
        id="manhattan_stacked empty list",
    ),
    pytest.param(
        lambda: _genomewide().plot_manhattan_qq_stacked([]),
        "GWAS DataFrame",
        id="manhattan_qq_stacked empty list",
    ),
    pytest.param(
        lambda: pylocuszoom.ManhattanPlotter(species=None).plot_manhattan(_GWAS),
        "custom_chrom_order",
        id="manhattan without chromosome order",
    ),
    pytest.param(
        lambda: pylocuszoom.ColocPlotter().plot_coloc(
            _COLOC, _COLOC.assign(pos=[1, 2])
        ),
        "overlapping",
        id="coloc without overlap",
    ),
    pytest.param(
        lambda: pylocuszoom.ColocPlotter().plot_coloc(_COLOC, _COLOC, ld_col="r2"),
        "r2",
        id="coloc missing ld column",
    ),
    pytest.param(
        lambda: pylocuszoom.ColocPlotter().plot_coloc(_COLOC, _COLOC, lead_snp="c"),
        "lead_snp",
        id="coloc lead_snp not found",
    ),
    pytest.param(
        lambda: pylocuszoom.LDHeatmapPlotter().plot_ld_heatmap(
            [[1.0, 0.5]], snp_ids=["a", "b"]
        ),
        "square",
        id="ld heatmap not square",
    ),
    pytest.param(
        lambda: pylocuszoom.LDHeatmapPlotter().plot_ld_heatmap(
            [[1.0, 0.5], [0.5, 1.0]], snp_ids=["a", "b"], lead_snp="c"
        ),
        "lead_snp",
        id="ld heatmap lead_snp not found",
    ),
    pytest.param(
        lambda: pylocuszoom.get_backend("nope"),
        "nope",
        id="unknown backend",
    ),
    pytest.param(
        lambda: pylocuszoom.load_gwas("results.txt", format="nope"),
        "nope",
        id="unknown loader format",
    ),
    pytest.param(
        lambda: _regional().plot([1, 2, 3], chrom=1, start=1, end=1000),
        "DataFrame",
        id="unsupported frame type",
    ),
]


class TestInputErrorsArePyLocusZoomErrors:
    """``except PyLocusZoomError`` catches every input error the public API raises."""

    @pytest.mark.parametrize(("call", "names"), INPUT_ERRORS)
    def test_input_error_is_a_validation_error_naming_the_input(self, call, names):
        with pytest.raises(ValidationError, match=names):
            call()
