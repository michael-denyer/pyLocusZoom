"""The table-driven engine every tabular loader runs through.

Most loaders share one shape: read a delimited file, rename source columns to
standard names, log a count, then validate. :class:`LoaderSpec` captures the
parts that differ between formats, and :func:`_load_tabular` is the single
engine that runs the shared steps. A new static format is one ``LoaderSpec``
constant plus a thin wrapper in the family module, not another copy-pasted
function body.
"""

import warnings
from collections.abc import Callable
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Any, Optional, Union

import pandas as pd

from ..logging import logger
from ..schemas import DEPRECATED_ALIAS_REMOVED_IN, Canonical
from ..utils import normalize_chrom_series
from ..validation import ColumnSpec, check

# Output-column tokens: a spec target equal to one of the `**out_cols` keys
# passed by the wrapper (`pos_col`, `p_col`, `rs_col`, `cs_col`) is substituted
# with the caller's chosen name; any other target is used literally. This is how
# the caller-configurable column names flow into an otherwise constant spec.


@dataclass(frozen=True)
class LoaderSpec:
    """Declarative description of how to load one tabular file format.

    Attributes:
        log_fmt: Debug log template; ``{n}`` is filled with the row count.
        read: Keyword arguments forwarded to ``pandas.read_csv``.
        schema: Validation to run after mapping. Takes the resolved output
            columns and returns the format's ``ColumnSpec``. None skips
            validation for a format that carries nothing checkable.
        schema_requires: Columns the format may structurally lack. An absent
            one is dropped from the schema's rules rather than failing, which
            is how CAVIAR and MatrixEQTL load without a position column.
        col_map: Static ``source -> target`` renames. Targets may be tokens.
        p_candidates: P-value source columns; the first present maps to ``p_col``.
        col_candidates: ``target -> candidates``; first present candidate maps to
            the (possibly token) target. Replaces hand-rolled first-match loops.
        transform: Optional ``(df, out_cols) -> df`` applied after renaming, for
            format-specific reshaping (e.g. credible-set assignment).
        gene_filter: ``"contains"`` or ``"exact"`` to enable eQTL gene filtering.
        clean_chrom: Strip a ``chr`` prefix from the ``chr`` column.
    """

    log_fmt: str
    read: dict[str, Any] = field(default_factory=dict)
    schema: Optional[Callable[[dict[str, str]], ColumnSpec]] = None
    schema_requires: tuple[str, ...] = ()
    col_map: dict[str, str] = field(default_factory=dict)
    p_candidates: tuple[str, ...] = ()
    col_candidates: dict[str, tuple[str, ...]] = field(default_factory=dict)
    transform: Optional[Callable[[pd.DataFrame, dict[str, str]], pd.DataFrame]] = None
    gene_filter: Optional[str] = None
    clean_chrom: bool = False


def _resolve(target: str, out_cols: dict[str, str]) -> str:
    """Substitute an output-column token, or return the literal target."""
    return out_cols.get(target, target)


def _first_present(df: pd.DataFrame, candidates: tuple[str, ...]) -> Optional[str]:
    """Return the first candidate column present in ``df``, else None."""
    for col in candidates:
        if col in df.columns:
            return col
    return None


def _filter_gene(df: pd.DataFrame, gene: str, mode: str) -> pd.DataFrame:
    """Filter eQTL rows to a gene by substring ("contains") or exact match."""
    if mode == "contains":
        return df[df["gene"].str.contains(gene, case=False, na=False)]
    return df[df["gene"] == gene]


def _relax(spec: ColumnSpec, absent: tuple[str, ...]) -> ColumnSpec:
    """Drop every rule naming a column the file structurally lacks."""
    if not absent:
        return spec
    return replace(
        spec,
        required=tuple(c for c in spec.required if c not in absent),
        numeric=tuple(c for c in spec.numeric if c not in absent),
        not_null=tuple(c for c in spec.not_null if c not in absent),
        ranges=tuple(r for r in spec.ranges if r.column not in absent),
    )


def _validate(df: pd.DataFrame, spec: LoaderSpec, out_cols: dict[str, str]) -> None:
    """Validate ``df`` against the format's schema.

    Every format validates strictly, so a failed column mapping is reported
    at load time where the source columns are still known. A column named in
    ``schema_requires`` is exempt only when the file genuinely lacks it.
    """
    if spec.schema is None:
        return
    absent = tuple(c for c in spec.schema_requires if c not in df.columns)
    check(df, _relax(spec.schema(out_cols), absent))


# The canonical name each output-column token resolves to when the caller
# names none. Overriding one is deprecated: a loader's job is to produce the
# package's vocabulary, and a caller who wants another name can rename.
_CANONICAL_OUT_COLS = {
    "pos_col": Canonical.POS,
    "p_col": Canonical.P,
    "rs_col": Canonical.RS,
}


def _canonical_out_cols(out_cols: dict[str, Optional[str]]) -> dict[str, str]:
    """Resolve the output-column tokens onto ``Canonical``, warning on overrides."""
    overridden = sorted(
        k for k, v in out_cols.items() if v is not None and k in _CANONICAL_OUT_COLS
    )
    if overridden:
        warnings.warn(
            f"{', '.join(overridden)} on the GWAS loaders is deprecated. Loaders "
            f"emit the canonical columns '{Canonical.CHROM}', '{Canonical.POS}', "
            f"'{Canonical.P}' and '{Canonical.RS}'; rename after loading instead. "
            f"The parameters are removed in {DEPRECATED_ALIAS_REMOVED_IN}.",
            DeprecationWarning,
            stacklevel=4,
        )
    return {k: _CANONICAL_OUT_COLS[k] if v is None else v for k, v in out_cols.items()}


def _load_tabular(
    filepath: Union[str, Path],
    spec: LoaderSpec,
    *,
    gene: Optional[str] = None,
    **requested: Optional[str],
) -> pd.DataFrame:
    """Load a tabular file per ``spec``: read, map, rename, transform, validate."""
    out_cols = _canonical_out_cols(requested)
    df = pd.read_csv(filepath, **spec.read)

    col_map = {src: _resolve(dst, out_cols) for src, dst in spec.col_map.items()}

    match = _first_present(df, spec.p_candidates)
    if match is not None:
        col_map[match] = out_cols["p_col"]

    for target, candidates in spec.col_candidates.items():
        match = _first_present(df, candidates)
        if match is not None:
            col_map[match] = _resolve(target, out_cols)

    df = df.rename(columns=col_map)

    if spec.transform is not None:
        df = spec.transform(df, out_cols)

    if spec.gene_filter is not None and gene is not None and "gene" in df.columns:
        df = _filter_gene(df, gene, spec.gene_filter)

    if spec.clean_chrom and "chr" in df.columns:
        df["chr"] = normalize_chrom_series(df["chr"])

    logger.debug(spec.log_fmt.format(n=len(df)))
    _validate(df, spec, out_cols)
    return df
