"""Shared plot-data intake policy."""

from dataclasses import dataclass
from types import MappingProxyType
from typing import List, Literal, Mapping, Optional, Type

import numpy as np
import pandas as pd

from .exceptions import ValidationError
from .logging import logger

P_VALUE_FLOOR = 1e-300
# Upper bound of the p-value domain, shared with ColumnSpec.pvalue so strict
# validation and plot-time intake cannot drift apart.
P_VALUE_MAX = 1.0

PValueFamily = Literal[
    "regional", "genome-wide", "qq", "eqtl", "phewas", "coloc", "loader"
]


@dataclass(frozen=True)
class PValuePolicy:
    """How one family treats a p-value that is null, non-numeric or out of range.

    Attributes:
        allow_zero: An exact zero is a valid p-value, clipped to
            ``P_VALUE_FLOOR``, rather than outside the ``(0, 1]`` domain.
        reject_invalid: An invalid p-value raises instead of dropping its row.
    """

    allow_zero: bool
    reject_invalid: bool


# The one statement of each family's p-value policy (ADR-0010). The families
# that draw many variants drop a few bad rows with a warning; the ones that
# draw one row per result or merge two sources reject them, because a dropped
# row there is a missing result.
P_VALUE_POLICY: Mapping[PValueFamily, PValuePolicy] = MappingProxyType(
    {
        # Regional association panels.
        "regional": PValuePolicy(allow_zero=True, reject_invalid=False),
        # Manhattan (including categorical), stacked and Miami panels.
        "genome-wide": PValuePolicy(allow_zero=True, reject_invalid=False),
        "qq": PValuePolicy(allow_zero=False, reject_invalid=False),
        # The regional eQTL panel.
        "eqtl": PValuePolicy(allow_zero=False, reject_invalid=False),
        "phewas": PValuePolicy(allow_zero=False, reject_invalid=True),
        "coloc": PValuePolicy(allow_zero=False, reject_invalid=True),
        # Every loader, through ColumnSpec.pvalue.
        "loader": PValuePolicy(allow_zero=False, reject_invalid=True),
    }
)


def pvalue_faults(values: pd.Series, column: str, *, allow_zero: bool) -> List[str]:
    """Describe every p-value in ``values`` outside the policy's domain.

    Args:
        values: The p-value column.
        column: Its name, for the messages.
        allow_zero: Whether an exact zero is valid.

    Returns:
        One message per kind of fault (null, non-numeric, below or above the
        domain), empty when every value is valid.
    """
    numeric = pd.to_numeric(values, errors="coerce")
    below = numeric < 0 if allow_zero else numeric <= 0
    counts = (
        (int(values.isna().sum()), f"Column '{column}' has {{}} null values"),
        (
            int(numeric.isna().sum() - values.isna().sum()),
            f"Column '{column}' has {{}} non-numeric values",
        ),
        (
            int(below.sum()),
            f"Column '{column}': {{}} values {'<' if allow_zero else '<='} 0",
        ),
        (
            int((numeric > P_VALUE_MAX).sum()),
            f"Column '{column}': {{}} values > {P_VALUE_MAX:g}",
        ),
    )
    return [message.format(count) for count, message in counts if count]


def prepare_pvalue_data(
    df: pd.DataFrame,
    p_col: str,
    family: PValueFamily,
    *,
    out_col: str = "neglog10p",
    on_empty: Optional[str] = None,
    error_class: Type[ValidationError] = ValidationError,
) -> pd.DataFrame:
    """Return a copy with valid p-values and a finite ``-log10(p)`` column.

    Null, non-numeric and out-of-range p-values are dropped or rejected as
    ``family``'s row of :data:`P_VALUE_POLICY` says. Tiny valid values are
    clipped before taking ``-log10``.

    Args:
        df: Frame holding the p-value column.
        p_col: Name of the p-value column.
        family: The plot family, which selects the policy.
        out_col: Name of the transformed column to write.
        on_empty: Message to raise when no row survives filtering. None
            returns the empty frame instead, which the regional path relies on.
        error_class: Exception raised when the policy rejects a p-value.

    Returns:
        A filtered copy with numeric ``p_col`` and the transformed ``out_col``.

    Raises:
        ValidationError: If the policy rejects an invalid p-value, or if
            nothing survives and ``on_empty`` names a message.
    """
    policy = P_VALUE_POLICY[family]
    if policy.reject_invalid:
        faults = pvalue_faults(df[p_col], p_col, allow_zero=policy.allow_zero)
        if faults:
            raise error_class(
                f"Invalid p-values for a {family} plot:\n"
                + "\n".join(f"  - {fault}" for fault in faults)
            )

    result = df.copy()
    initial_count = len(result)
    p_values = pd.to_numeric(result[p_col], errors="coerce")
    lower_mask = p_values >= 0 if policy.allow_zero else p_values > 0
    valid = p_values.notna() & lower_mask & (p_values <= P_VALUE_MAX)
    dropped = int((~valid).sum())
    nan_count = int(result[p_col].isna().sum())
    if nan_count:
        logger.warning("Found {} NaN p-values, filtering out", nan_count)
    numeric_missing = int(p_values.isna().sum()) - nan_count
    if numeric_missing:
        logger.warning("Found {} non-numeric p-values, filtering out", numeric_missing)
    out_of_range = int((p_values.notna() & ~valid).sum())
    if out_of_range:
        logger.warning(
            "Found {} p-values outside {} range, filtering out",
            out_of_range,
            "[0, 1]" if policy.allow_zero else "(0, 1]",
        )
    result = result.loc[valid].copy()
    valid_values = p_values.loc[valid]
    result[p_col] = valid_values
    clipped = int((valid_values < P_VALUE_FLOOR).sum())
    if clipped:
        logger.debug("Clipping {} p-values below {}", clipped, P_VALUE_FLOOR)
    result[out_col] = -np.log10(valid_values.clip(lower=P_VALUE_FLOOR))
    if dropped:
        logger.debug("P-value filtering removed {} of {} rows", dropped, initial_count)
    if result.empty and on_empty is not None:
        raise ValidationError(on_empty)
    return result
