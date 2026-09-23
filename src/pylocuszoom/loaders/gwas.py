"""GWAS summary-statistic loaders, and the front door that detects the format.

Every ``load_*`` GWAS loader takes the same arguments and returns the same
shape, so the per-loader docstrings below cover only what is specific to their
format:

Args:
    filepath: Path to the results file.

Returns:
    DataFrame in the canonical column vocabulary: chr, pos, p_value, rs.
    Rename columns after loading if you want other names.

Raises:
    LoaderValidationError: If the format's columns cannot be mapped or the
        mapped values fail the format's schema.
"""

from collections.abc import Callable
from pathlib import Path
from typing import Optional, Union

import pandas as pd

from .._data import P_VALUE_FLOOR
from ..exceptions import ValidationError
from ..logging import logger
from ..schemas import Canonical, gwas_load_spec
from ._engine import LoaderSpec, _load_tabular

_GWAS_SCHEMA = gwas_load_spec()


# No comment= here on purpose. PLINK 2's --glm header line starts with "#CHROM",
# so a "#" comment prefix makes pandas discard the header and promote the first
# variant to column labels. PLINK writes no comment lines of its own.
_PLINK_SPEC = LoaderSpec(
    log_fmt="Loaded PLINK file with {n} variants",
    read={"sep": r"\s+"},
    col_candidates={
        Canonical.POS: ("BP", "POS", "bp", "pos"),
        Canonical.P: ("P", "P_BOLT_LMM", "p", "PVAL", "pval", "P_LINREG"),
        Canonical.RS: ("SNP", "ID", "rsid", "RSID", "MarkerName", "variant_id"),
        Canonical.CHROM: ("CHR", "chr", "CHROM", "chrom", "#CHROM"),
    },
    schema=_GWAS_SCHEMA,
)


def load_plink_assoc(filepath: Union[str, Path]) -> pd.DataFrame:
    """Load PLINK association results (.assoc, .assoc.linear, .assoc.logistic, .qassoc).

    Automatically detects PLINK format variant and maps columns to standard names.
    See the module docstring for the shared GWAS loader arguments and return value.

    Example:
        >>> gwas_df = load_plink_assoc("results.assoc.linear")
        >>> fig = plotter.plot(gwas_df, chrom=1, start=1e6, end=2e6)
    """
    return _load_tabular(filepath, _PLINK_SPEC)


def _regenie_pvalue(df: pd.DataFrame, out_cols: dict[str, str]) -> pd.DataFrame:
    """Derive REGENIE p-values: prefer computed LOG10P, else rename P."""
    p_col = Canonical.P
    if "LOG10P" in df.columns:
        # REGENIE writes LOG10P so that p below the float range survives. The
        # plot clips p at P_VALUE_FLOOR anyway; clipping here keeps those rows
        # positive instead of underflowing to 0 and failing the p > 0 check.
        df[p_col] = (10.0 ** -df["LOG10P"].astype(float)).clip(lower=P_VALUE_FLOOR)
    elif "P" in df.columns:
        df = df.rename(columns={"P": p_col})
    return df


_REGENIE_SPEC = LoaderSpec(
    log_fmt="Loaded REGENIE file with {n} variants",
    read={"sep": r"\s+", "comment": "#"},
    col_map={"GENPOS": Canonical.POS, "ID": Canonical.RS, "CHROM": Canonical.CHROM},
    transform=_regenie_pvalue,
    schema=_GWAS_SCHEMA,
)


def load_regenie(filepath: Union[str, Path]) -> pd.DataFrame:
    """Load REGENIE association results (.regenie).

    Prefers the computed LOG10P column over a raw P column.
    See the module docstring for the shared GWAS loader arguments and return value.

    Example:
        >>> gwas_df = load_regenie("results.regenie")
    """
    return _load_tabular(filepath, _REGENIE_SPEC)


_BOLT_LMM_SPEC = LoaderSpec(
    log_fmt="Loaded BOLT-LMM file with {n} variants",
    read={"sep": "\t"},
    col_map={"BP": Canonical.POS, "SNP": Canonical.RS, "CHR": Canonical.CHROM},
    col_candidates={Canonical.P: ("P_BOLT_LMM", "P_BOLT_LMM_INF")},
    schema=_GWAS_SCHEMA,
)


def load_bolt_lmm(filepath: Union[str, Path]) -> pd.DataFrame:
    """Load BOLT-LMM association results (.stats).

    See the module docstring for the shared GWAS loader arguments and return value.

    Example:
        >>> gwas_df = load_bolt_lmm("results.stats")
    """
    return _load_tabular(filepath, _BOLT_LMM_SPEC)


_GEMMA_SPEC = LoaderSpec(
    log_fmt="Loaded GEMMA file with {n} variants",
    read={"sep": "\t"},
    col_map={"ps": Canonical.POS, "rs": Canonical.RS, "chr": Canonical.CHROM},
    col_candidates={Canonical.P: ("p_wald", "p_lrt", "p_score")},
    schema=_GWAS_SCHEMA,
)


def load_gemma(filepath: Union[str, Path]) -> pd.DataFrame:
    """Load GEMMA association results (.assoc.txt).

    See the module docstring for the shared GWAS loader arguments and return value.

    Example:
        >>> gwas_df = load_gemma("output.assoc.txt")
    """
    return _load_tabular(filepath, _GEMMA_SPEC)


_SAIGE_SPEC = LoaderSpec(
    log_fmt="Loaded SAIGE file with {n} variants",
    read={"sep": "\t"},
    col_map={"POS": Canonical.POS, "MarkerID": Canonical.RS, "CHR": Canonical.CHROM},
    col_candidates={Canonical.P: ("p.value.NA", "p.value")},
    schema=_GWAS_SCHEMA,
)


def load_saige(filepath: Union[str, Path]) -> pd.DataFrame:
    """Load SAIGE association results.

    See the module docstring for the shared GWAS loader arguments and return value.

    Example:
        >>> gwas_df = load_saige("results.txt")
    """
    return _load_tabular(filepath, _SAIGE_SPEC)


_GWAS_CATALOG_SPEC = LoaderSpec(
    log_fmt="Loaded GWAS Catalog file with {n} variants",
    read={"sep": "\t"},
    col_map={
        "base_pair_location": Canonical.POS,
        "variant_id": Canonical.RS,
        "chromosome": Canonical.CHROM,
        "p_value": Canonical.P,
    },
    schema=_GWAS_SCHEMA,
)


def load_gwas_catalog(filepath: Union[str, Path]) -> pd.DataFrame:
    """Load GWAS Catalog summary statistics format.

    See the module docstring for the shared GWAS loader arguments and return value.
    """
    return _load_tabular(filepath, _GWAS_CATALOG_SPEC)


# Filename substrings that identify a GWAS format.
_FORMAT_HINTS: tuple[tuple[str, str], ...] = (
    (".assoc", "plink"),
    (".qassoc", "plink"),
    (".glm.", "plink"),
    (".regenie", "regenie"),
    (".stats", "bolt"),
    ("gemma", "gemma"),
    (".assoc.txt", "gemma"),
    ("saige", "saige"),
)

_GWAS_LOADERS: dict[str, Callable[..., pd.DataFrame]] = {
    "plink": load_plink_assoc,
    "regenie": load_regenie,
    "bolt": load_bolt_lmm,
    "gemma": load_gemma,
    "saige": load_saige,
    "catalog": load_gwas_catalog,
}


def _detect_format(filepath: Path) -> str:
    """Identify a GWAS format from the filename, defaulting to PLINK."""
    name = filepath.name.lower()
    for hint, format in sorted(_FORMAT_HINTS, key=lambda h: len(h[0]), reverse=True):
        if hint in name:
            return format
    logger.warning(
        f"Could not detect GWAS file format for '{filepath}'. "
        "Defaulting to 'plink'. Specify format= parameter explicitly."
    )
    return "plink"


def load_gwas(
    filepath: Union[str, Path],
    format: Optional[str] = None,
) -> pd.DataFrame:
    """Load GWAS results with automatic format detection.

    Args:
        filepath: Path to GWAS results file.
        format: File format. If None, auto-detects from extension.
            Options: "plink", "regenie", "bolt", "gemma", "saige", "catalog".

    Returns:
        DataFrame in the canonical column vocabulary: chr, pos, p_value, rs.

    Raises:
        ValidationError: If ``format`` names an unknown format.

    Example:
        >>> # Auto-detect format
        >>> gwas_df = load_gwas("results.assoc.linear")
        >>>
        >>> # Explicit format
        >>> gwas_df = load_gwas("results.txt", format="regenie")
    """
    filepath = Path(filepath)
    format = format or _detect_format(filepath)

    if format not in _GWAS_LOADERS:
        raise ValidationError(
            f"Unknown format '{format}'. Options: {list(_GWAS_LOADERS.keys())}"
        )

    return _GWAS_LOADERS[format](filepath)
