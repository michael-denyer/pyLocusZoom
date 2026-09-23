"""DataFrame validation engine for pyLocusZoom.

Declares the rule vocabulary (:class:`ColumnSpec`, :class:`RangeRule`) and
runs it against a frame with :func:`check`, accumulating every fault before
raising. The families the rules describe live in ``schemas.py``.
"""

import operator
from dataclasses import dataclass
from typing import List, Optional, Set, Tuple, Type

import pandas as pd
from pandas.api.types import is_numeric_dtype

from ._data import P_VALUE_POLICY, pvalue_faults
from .exceptions import ValidationError

_COMPARE = {
    "<": operator.lt,
    "<=": operator.le,
    ">": operator.gt,
    ">=": operator.ge,
}


@dataclass(frozen=True)
class RangeRule:
    """A single numeric-range constraint on one column."""

    column: str
    min_val: Optional[float] = None
    max_val: Optional[float] = None
    exclusive_min: bool = False
    exclusive_max: bool = False


@dataclass(frozen=True)
class ColumnSpec:
    """Declarative validation contract for one DataFrame family.

    Each field names the columns a rule applies to. :func:`check` runs the
    rules in a fixed order and raises ``error_class`` if any fails.

    Args:
        name: Dataset name used in error messages (e.g. "GWAS").
        required: Columns that must exist.
        numeric: Columns that must have a numeric dtype.
        not_null: Columns that must contain no nulls.
        ranges: Numeric-range constraints applied in order.
        pvalue: The p-value column. It must be numeric, and every value
            null, non-numeric or outside the domain is rejected, as the
            ``"loader"`` row of ``_data.P_VALUE_POLICY`` says; do not list it
            under ``numeric`` or ``not_null`` as well.
        ordering: ``(lower, upper)`` pairs where lower must never exceed upper.
        non_empty: Reject a frame with no rows, before any column rule runs.
        error_class: Exception raised on failure.
    """

    name: str
    required: Tuple[str, ...] = ()
    numeric: Tuple[str, ...] = ()
    not_null: Tuple[str, ...] = ()
    ranges: Tuple[RangeRule, ...] = ()
    pvalue: Optional[str] = None
    ordering: Tuple[Tuple[str, str], ...] = ()
    non_empty: bool = False
    error_class: Type[ValidationError] = ValidationError


def _range_errors(
    df: pd.DataFrame, rule: RangeRule, non_numeric: Set[str]
) -> List[str]:
    """Report the range faults in one column, skipping what cannot be compared.

    pandas raises ``TypeError`` comparing a non-numeric column to a numeric
    bound, which would mask the structured error.
    """
    if rule.column not in df.columns or rule.column in non_numeric:
        return []

    values = df[rule.column].dropna()
    bounds = (
        (rule.min_val, "<=" if rule.exclusive_min else "<"),
        (rule.max_val, ">=" if rule.exclusive_max else ">"),
    )
    return [
        f"Column '{rule.column}': {count} values {symbol} {bound}"
        for bound, symbol in bounds
        if bound is not None and (count := _COMPARE[symbol](values, bound).sum()) > 0
    ]


def resolve_column(
    df: pd.DataFrame,
    column: Optional[str],
    *,
    parameter: str,
    optional_default: Optional[str] = None,
    frame: str = "the data",
    error_class: Type[ValidationError] = ValidationError,
) -> Optional[str]:
    """Resolve a column the caller named, or None when the feature is off.

    The one owner of "named means required". A column the caller names must
    be in the frame. The exception is a canonical default for an optional
    feature (``rs``, ``cs``, ``category``): left at that default, a frame
    without the column simply draws without the feature.

    Args:
        df: The frame the column should be in.
        column: The caller's column name, or None for the feature off.
        parameter: The option that named the column, for the error message.
        optional_default: The canonical name that may be absent, if any.
        frame: What the frame is, for the error message.
        error_class: Exception raised for a missing named column.

    Returns:
        ``column`` if the frame carries it, else None when ``column`` is None
        or the absent optional default.

    Raises:
        ValidationError: If ``column`` names a column the frame lacks.
    """
    if column is None or column in df.columns:
        return column
    if column == optional_default:
        return None
    raise error_class(
        f"{parameter}='{column}' names no column of {frame}. "
        f"Available: {list(df.columns)}"
    )


def check(df: pd.DataFrame, spec: ColumnSpec) -> None:
    """Validate ``df`` against ``spec``, accumulating all faults before raising.

    Args:
        df: DataFrame to validate.
        spec: The validation contract to apply.

    Raises:
        ValidationError: If any rule fails. The concrete type is
            ``spec.error_class``.
    """
    errors: List[str] = []
    non_numeric: Set[str] = set()

    if spec.non_empty and df.empty:
        errors.append("is empty — no rows to plot")

    missing = [col for col in spec.required if col not in df.columns]
    if missing:
        errors.append(f"Missing columns: {missing}. Available: {list(df.columns)}")

    numeric = spec.numeric if spec.pvalue is None else (*spec.numeric, spec.pvalue)
    for col in numeric:
        if col in df.columns and not is_numeric_dtype(df[col]):
            errors.append(f"Column '{col}' must be numeric, got {df[col].dtype}")
            non_numeric.add(col)

    for col in spec.not_null:
        if col in df.columns:
            null_count = df[col].isna().sum()
            if null_count > 0:
                errors.append(f"Column '{col}' has {null_count} null values")

    for rule in spec.ranges:
        errors.extend(_range_errors(df, rule, non_numeric))

    if spec.pvalue in df.columns and spec.pvalue not in non_numeric:
        errors.extend(
            pvalue_faults(
                df[spec.pvalue],
                spec.pvalue,
                allow_zero=P_VALUE_POLICY["loader"].allow_zero,
            )
        )

    for lower_col, upper_col in spec.ordering:
        if any(
            col not in df.columns or col in non_numeric
            for col in (lower_col, upper_col)
        ):
            continue
        inverted = (df[lower_col] > df[upper_col]).sum()
        if inverted > 0:
            errors.append(f"{inverted} rows have {lower_col} > {upper_col}")

    if errors:
        raise spec.error_class(
            f"{spec.name} validation failed:\n"
            + "\n".join(f"  - {error}" for error in errors)
        )
